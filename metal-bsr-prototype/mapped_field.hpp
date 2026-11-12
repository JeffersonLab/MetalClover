#pragma once
#include <fcntl.h>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace MetalCloverEnv {

  // File-backed FP64 storage. Explicit slab writeback/drop bounds resident
  // pages; virtual addresses stay valid for halo reads in the full Wilson
  // stencil.
  struct MappedField {
    int fd = -1;
    size_t bytes = 0;
    void *data = nullptr;
    bool writable = false;
    std::string path;
    MappedField(const std::string &name, size_t size = 0, bool create = false)
        : path(name) {
      if (name.empty()) {
        bytes = size;
        writable = true;
        data = mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                    MAP_PRIVATE | MAP_ANON, -1, 0);
        if (data == MAP_FAILED)
          throw std::runtime_error("anonymous mapping");
        return;
      }
      fd = ::open(name.c_str(),
                  create ? O_RDWR | O_CREAT | O_EXCL
                         : (size ? O_RDWR : O_RDONLY),
                  0600);
      if (fd < 0)
        throw std::runtime_error("open mapped field: " + name);
      writable = create || size;
      struct stat s{};
      if (fstat(fd, &s))
        throw std::runtime_error("stat mapped field");
      if (create && ftruncate(fd, size))
        throw std::runtime_error("size mapped field");
      bytes = create ? size : size_t(s.st_size);
      if (size && !create && bytes != size)
        throw std::runtime_error("mapped field size mismatch");
      data = mmap(nullptr, bytes, PROT_READ | (writable ? PROT_WRITE : 0),
                  MAP_SHARED, fd, 0);
      if (data == MAP_FAILED)
        throw std::runtime_error("mmap failed");
      madvise(data, bytes, MADV_RANDOM);
    }
    MappedField(const MappedField &) = delete;
    ~MappedField() {
      if (data && data != MAP_FAILED)
        munmap(data, bytes);
      if (fd >= 0)
        close(fd);
    }
    void drop(size_t offset = 0, size_t length = 0, bool writeback = false) {
      if (!length)
        length = bytes - offset;
      size_t page = size_t(sysconf(_SC_PAGESIZE)), first = offset / page * page,
             last = std::min(bytes, (offset + length + page - 1) / page * page);
      if (writeback && fd >= 0 &&
          msync(static_cast<char *>(data) + first, last - first, MS_SYNC))
        throw std::runtime_error("FP64 slab writeback failed");
      madvise(static_cast<char *>(data) + first, last - first,
              fd >= 0 ? MADV_DONTNEED : MADV_FREE);
    }
    void rename_to(const std::string &destination) {
      if (std::filesystem::exists(destination))
        throw std::runtime_error("output exists");
      if (::rename(path.c_str(), destination.c_str()))
        throw std::runtime_error("rename solved field failed");
      path = destination;
    }
  };

} // namespace MetalCloverEnv
