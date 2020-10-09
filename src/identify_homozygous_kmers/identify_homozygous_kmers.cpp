#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>

#include <jellyfish/rectangular_binary_matrix.hpp>
#include <jellyfish/cpp_array.hpp>
#include <jellyfish/err.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/mer_heap.hpp>
#include <jellyfish/misc.hpp>

namespace err = jellyfish::err;

using jellyfish::file_header;
using jellyfish::RectangularBinaryMatrix;
using jellyfish::mer_dna;
using jellyfish::cpp_array;
typedef std::unique_ptr<binary_reader> binary_reader_ptr;
typedef std::unique_ptr<text_reader> text_reader_ptr;

struct file_info {
    std::ifstream is;
    file_header header;

    file_info(const char* path):
        is(path),
        header(is)
    { }
};

struct common_info {
    unsigned int key_len;
    size_t max_reprobe_offset;
    size_t size;
    unsigned int out_counter_len;
    std::string format;
    RectangularBinaryMatrix matrix;

    common_info(RectangularBinaryMatrix&& m) : matrix(std::move(m))
    { }
};

common_info read_headers(char* read_file, char* assembly_file, cpp_array<file_info>& files) {
    files.init(0, read_file);
    if (!files[0].is.good()) {
        err::die(err::msg() << "failed to open input file '" << read_file << "'");
    }

    file_header& h = files[0].header;
    common_info res(h.matrix());
    res.key_len = h.key_len();
    res.max_reprobe_offset = h.max_reprobe_offset();
    res.size = h.size();
    res.format = h.format();
    res.out_counter_len = h.counter_len();

    files.init(1, assembly_file);
    file_header& nh = files[1].header;
    if (!files[1].is.good()) {
        err::die(err::msg() << "failed to open input file '" << assembly_file << "'");
    }
    if (res.format != nh.format()) {
        err::die(err::msg() << "can't compare files with different formats (" << res.format << ", " << nh.format() << ")");
    }
    if (res.key_len != nh.key_len()) {
        err::die(err::msg() << "can't compare hashes of different key lengths (" << res.key_len << ", " << nh.key_len() << ")");
    }
    if (res.max_reprobe_offset != nh.max_reprobe_offset()) {
        err::die(err::msg() << "can't compare hashes with different reprobing strategies");
    }
    if (res.size != nh.size()) {
        err::die(err::msg() << "can't compare hash with different size (" << res.size << ", " << nh.size() << ")");
    }
    if (res.matrix != nh.matrix()) {
        err::die(err::msg() << "can't compare hash with different hash function");
    }

    return res;
}

template<typename reader_type>
void output_counts(cpp_array<file_info>& files, int low_read_count, int high_read_count) {
    cpp_array<reader_type> readers(files.size());
    typedef jellyfish::mer_heap::heap<mer_dna, reader_type> heap_type;
    typedef typename heap_type::const_item_t heap_item;
    heap_type heap(files.size());

    for (size_t i = 0; i < files.size(); ++i) {
        readers.init(i, files[i].is, &files[i].header);
        if (readers[i].next()) {
            heap.push(readers[i]);
        }
    }

    heap_item head = heap.head();
    mer_dna key;
    const int num_files = files.size();
    const reader_type* base = &readers[0];
    uint64_t counts[num_files];
    int kmer_count = 0;

    while (heap.is_not_empty()) {
        key = head->key_;
        memset(counts, '\0', sizeof(uint64_t) * num_files);
        do {
            counts[head->it_ - base] = head->val_;
            heap.pop();
            if (head->it_->next()) {
                heap.push(*head->it_);
            }
            head = heap.head();
        } while (head->key_ == key && heap.is_not_empty());

        if (counts[0] <= high_read_count && counts[0] >= low_read_count && counts[1] == 2) {
            ++kmer_count;
            // std::cout << key << " " << counts[0] << " " << counts[1] << "\n";
            std::cout << ">kmer" << kmer_count << "\n" << key << "\n";
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc < 4) {err::die(err::msg() << "usage: " <<
        argv[0] << "hom_peak_low hom_peak_high read_file assembly_file" <<
        "\n\nARGUMENTS:\n" <<
        "\thom_peak_low\t\tk-mer coverage in read data that represents the lower end of the homozygous peak\n" <<
        "\thom_peak_high\t\tk-mer coverage in read data that represents the high end of the homozygous peak\n" <<
        "\tread_file\t\tjellyfish database from short read data\n" <<
        "\tassembly_file\t\tjellyfish database from genome assembly\n");
    }

    unsigned int hom_peak_low = atoi(argv[1]);
    unsigned int hom_peak_high = atoi(argv[2]);

    cpp_array<file_info> files(2);
    common_info cinfo = read_headers(argv[3], argv[4], files);
    mer_dna::k(cinfo.key_len / 2);

    if (cinfo.format == binary_dumper::format) {
        output_counts<binary_reader>(files, hom_peak_low, hom_peak_high);
    } else if (cinfo.format == text_dumper::format) {
        output_counts<text_reader>(files, hom_peak_low, hom_peak_high);
    } else {
        err::die(err::msg() << "format '" << cinfo.format << "' not supported\n");
    }

    return 0;
}
