#include "Filesystem.hpp"

#include <string>
#include <fstream>
#include <stdio.h>  // defines FILENAME_MAX
#include <cctype>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>    // errno, ENOENT, EEXIST

#ifdef _WIN32
  #include <direct.h>   // _mkdir
  #define getCurrentDir _getcwd
#else
  #include <unistd.h>
  #define getCurrentDir getcwd
#endif

#include "Output.hpp"

namespace Shtns { namespace filesystem {

std::string path::string() const
{
  if (empty())
    return "";

  auto it = begin();
  auto result = *it;
  for (++it; it != end(); ++it)
    result += preferred_separator + *it;
  return result;
}


void path::split(std::string p)
{
  std::string separators = "/\\";
  bool relative = true;

  auto pos = p.find_first_of(separators);
  if (is_absolute(p))
    pos = p.find_first_of(separators,pos+1);
  for (; pos != std::string::npos;
         pos = p.find_first_of(separators))
  {
    auto token = p.substr(0,pos);
    trim(token);
    if (!token.empty() && token != "." && token != "..") {
      this->push_back(token);
      relative = false;
    } else if (token == "..") {
      if (relative || this->empty())
        this->push_back(token);
      else
        this->pop_back();
    }
    p = p.substr(pos+1);
  }
  this->push_back(p);
}


path path::stem() const
{
  auto f = filename().string();
  auto pos = f.find_last_of('.');
  if (f == "." || f == ".." || pos == std::string::npos)
    return {f};
  else
    return {f.substr(0,pos)};
}


path path::extension() const
{
  auto f = filename().string();
  auto pos = f.find_last_of('.');
  if (f == "." || f == ".." || pos == std::string::npos)
    return {};
  else
    return {f.substr(pos)};
}


bool path::is_absolute(std::string p)
{
  if (p[0] == '/')
    return true;

  // c:\ or z:/
  if (std::isalpha(p[0]) && p[1] == ':' && (p[2] == '/' || p[2] == '\\'))
    return true;

  return false;
}


path& path::operator/=(path const& p)
{
  insert(end(), p.begin(), p.end());
  original += preferred_separator + p.original;
  return *this;
}


bool path::is_file() const
{
  std::string p = this->string();
  struct stat info;
  return stat(p.c_str(), &info) == 0 && (info.st_mode & S_IFREG) != 0;
}

bool path::is_directory() const
{
  std::string p = this->string();
  struct stat info;
  return stat(p.c_str(), &info) == 0 && (info.st_mode & S_IFDIR) != 0;
}


path current_path()
{
  char cwd_[FILENAME_MAX];
  [[maybe_unused]] auto _ =  getCurrentDir(cwd_, sizeof(cwd_));
  std::string cwd(cwd_);
  return { trim(cwd) };
}


bool exists(path const& p)
{
  return p.is_file() || p.is_directory();
}


bool create_directories(path const& p)
{
  if (p.is_directory())
    return true;

  auto parent = p.parent_path();
  if (!parent.empty() && !parent.is_directory())
    create_directories(parent);

#ifdef _WIN32
  int ret = _mkdir(p.string().c_str());
#else
  mode_t mode = 0755;
  int ret = mkdir(p.string().c_str(), mode);
#endif
  if (ret == 0)
    return true;

  switch (errno)
  {
    case ENOENT:
      error_exit("parent didn't exist. Should not happen, since parent directory created before!");
      return false;
      break;
    case EEXIST:
      return true;
      break;
    default:
      return false;
  }
}

} } // end namespace Shtns
