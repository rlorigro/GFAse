#ifndef GFASE_BINARYIO_CPP_H
#define GFASE_BINARYIO_CPP_H

#include "BinaryIO.hpp"


void write_string_to_binary(ostream& file, string& s){
    ///
    /// Without worrying about size conversions, write any string to a file using ostream.write
    ///

    file.write(reinterpret_cast<const char*>(s.data()), s.size());
}


void read_string_from_binary(istream& file, string& s, uint64_t length){
    ///
    /// Without worrying about size conversions, read any value from a file using istream.read
    ///

    s.resize(length);
    file.read(reinterpret_cast<char*>(s.data()), length);
}


void pread_bytes(int file_descriptor, char* buffer_pointer, size_t bytes_to_read, off_t& offset){
    ///
    /// Reimplementation of binary read_bytes(), but with Linux pread, which is threadsafe
    ///

    while (bytes_to_read) {
        const ssize_t byte_count = ::pread(file_descriptor, buffer_pointer, bytes_to_read, offset);
        if (byte_count <= 0) {
            throw runtime_error("ERROR " + std::to_string(errno) + " while reading: " + string(::strerror(errno)));
        }
        bytes_to_read -= byte_count;
        buffer_pointer += byte_count;
        offset += byte_count;
    }
}


void pread_string_from_binary(int file_descriptor, string& s, uint64_t length, off_t& offset){
    ///
    /// Reimplementation of binary read_string_from_binary(), but with Linux pread, which is threadsafe
    ///

    s.resize(length);

    size_t bytes_to_read = length;
    char* buffer_pointer = reinterpret_cast<char*>(s.data());

    pread_bytes(file_descriptor, buffer_pointer, bytes_to_read, offset);
}

#endif //GFASE_BINARYIO_CPP_H
