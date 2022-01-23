#include "io.h"

#include <fcntl.h>
#if defined(__APPLE__)
#else
#include <malloc.h>
#endif
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <fstream>

namespace gbbs {
namespace gbbs_io {

#ifdef PARLAY_USE_STD_ALLOC
__mallopt::__mallopt() {
  mallopt(M_MMAP_MAX, 0);
  mallopt(M_TRIM_THRESHOLD, -1);
}

__mallopt __mallopt_var = __mallopt();
#endif

// returns a pointer and a length
std::pair<char*, size_t> mmapStringFromFile(const char* filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG(sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char* p =
      static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
  return std::make_pair(p, n);
}

void unmmap(const char* bytes, size_t bytes_size) {
  if (bytes) {
    const void* b = bytes;
    if (munmap(const_cast<void*>(b), bytes_size) == -1) {
      perror("munmap");
      exit(-1);
    }
  }
}

sequence<char> readStringFromFile(const char* fileName) {
  std::ifstream file(fileName, std::ios::in | std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    debug(std::cout << "# Unable to open file: " << fileName << "\n";);
    abort();
  }
  uint64_t end = file.tellg();
  file.seekg(0, std::ios::beg);
  uint64_t n = end - file.tellg();
  auto bytes = sequence<char>(n);  // n+1?
  file.read(bytes.begin(), n);
  file.close();
  return bytes;
}

std::tuple<char*, size_t> read_o_direct(const char* fname) {
  /* read using O_DIRECT, which bypasses caches. */
  int fd;
#if defined(__APPLE__)
  if ((fd = open(fname, O_RDONLY)) != -1) {
#else
  if ((fd = open(fname, O_RDONLY | O_DIRECT)) != -1) {
#endif
    debug(std::cout << "# input opened!"
                    << "\n";);
  } else {
    std::cout << "# can't open input file!";
  }
  //    posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED);

  size_t fsize = lseek(fd, 0, SEEK_END);
  lseek(fd, 0, 0);

/* allocate properly memaligned buffer for bytes */
#if defined(__APPLE__)
  char* bytes = NULL;
  posix_memalign((void**)&bytes, 4096 * 2, fsize + 4096);
#else
  char* bytes = (char*)memalign(4096 * 2, fsize + 4096);
#endif
  debug(std::cout << "# fsize = " << fsize << "\n";);

  size_t sz = 0;

  size_t pgsize = getpagesize();
  debug(std::cout << "# pgsize = " << pgsize << "\n";);

  size_t read_size = 1024 * 1024 * 1024;
  if (sz + read_size > fsize) {
    size_t k = std::ceil((fsize - sz) / pgsize);
    read_size = std::max(k * pgsize, pgsize);
    debug(std::cout << "# set read size to: " << read_size << " "
                    << (fsize - sz) << " bytes left"
                    << "\n";);
  }

  while (sz + read_size < fsize) {
    void* buf = bytes + sz;
    debug(std::cout << "# reading: " << read_size << "\n";);
    sz += read(fd, buf, read_size);
    debug(std::cout << "# read: " << sz << " bytes"
                    << "\n";);
    if (sz + read_size > fsize) {
      size_t k = std::ceil((fsize - sz) / pgsize);
      read_size = std::max(k * pgsize, pgsize);
      debug(std::cout << "# set read size to: " << read_size << " "
                      << (fsize - sz) << " bytes left"
                      << "\n";);
    }
  }
  if (sz < fsize) {
    debug(std::cout << "# last read: rem = " << (fsize - sz) << "\n";);
    void* buf = bytes + sz;
    sz += read(fd, buf, pgsize);
    debug(std::cout << "# read " << sz << " bytes "
                    << "\n";);
  }
  close(fd);
  return std::make_tuple(bytes, fsize);
}

}  // namespace gbbs_io
}  // namespace gbbs
