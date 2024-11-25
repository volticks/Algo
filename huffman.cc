#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <map>
#include <memory>
#include <unordered_map>
#include <vector>

// https://stackoverflow.com/questions/2602823/in-c-c-whats-the-simplest-way-to-reverse-the-order-of-bits-in-a-byte
// Cheesy, i know
unsigned char reverse_byte(unsigned char x, uint8_t n) {
  static const unsigned char table[] = {
      0x00, 0x80, 0x40, 0xc0, 0x20, 0xa0, 0x60, 0xe0, 0x10, 0x90, 0x50, 0xd0,
      0x30, 0xb0, 0x70, 0xf0, 0x08, 0x88, 0x48, 0xc8, 0x28, 0xa8, 0x68, 0xe8,
      0x18, 0x98, 0x58, 0xd8, 0x38, 0xb8, 0x78, 0xf8, 0x04, 0x84, 0x44, 0xc4,
      0x24, 0xa4, 0x64, 0xe4, 0x14, 0x94, 0x54, 0xd4, 0x34, 0xb4, 0x74, 0xf4,
      0x0c, 0x8c, 0x4c, 0xcc, 0x2c, 0xac, 0x6c, 0xec, 0x1c, 0x9c, 0x5c, 0xdc,
      0x3c, 0xbc, 0x7c, 0xfc, 0x02, 0x82, 0x42, 0xc2, 0x22, 0xa2, 0x62, 0xe2,
      0x12, 0x92, 0x52, 0xd2, 0x32, 0xb2, 0x72, 0xf2, 0x0a, 0x8a, 0x4a, 0xca,
      0x2a, 0xaa, 0x6a, 0xea, 0x1a, 0x9a, 0x5a, 0xda, 0x3a, 0xba, 0x7a, 0xfa,
      0x06, 0x86, 0x46, 0xc6, 0x26, 0xa6, 0x66, 0xe6, 0x16, 0x96, 0x56, 0xd6,
      0x36, 0xb6, 0x76, 0xf6, 0x0e, 0x8e, 0x4e, 0xce, 0x2e, 0xae, 0x6e, 0xee,
      0x1e, 0x9e, 0x5e, 0xde, 0x3e, 0xbe, 0x7e, 0xfe, 0x01, 0x81, 0x41, 0xc1,
      0x21, 0xa1, 0x61, 0xe1, 0x11, 0x91, 0x51, 0xd1, 0x31, 0xb1, 0x71, 0xf1,
      0x09, 0x89, 0x49, 0xc9, 0x29, 0xa9, 0x69, 0xe9, 0x19, 0x99, 0x59, 0xd9,
      0x39, 0xb9, 0x79, 0xf9, 0x05, 0x85, 0x45, 0xc5, 0x25, 0xa5, 0x65, 0xe5,
      0x15, 0x95, 0x55, 0xd5, 0x35, 0xb5, 0x75, 0xf5, 0x0d, 0x8d, 0x4d, 0xcd,
      0x2d, 0xad, 0x6d, 0xed, 0x1d, 0x9d, 0x5d, 0xdd, 0x3d, 0xbd, 0x7d, 0xfd,
      0x03, 0x83, 0x43, 0xc3, 0x23, 0xa3, 0x63, 0xe3, 0x13, 0x93, 0x53, 0xd3,
      0x33, 0xb3, 0x73, 0xf3, 0x0b, 0x8b, 0x4b, 0xcb, 0x2b, 0xab, 0x6b, 0xeb,
      0x1b, 0x9b, 0x5b, 0xdb, 0x3b, 0xbb, 0x7b, 0xfb, 0x07, 0x87, 0x47, 0xc7,
      0x27, 0xa7, 0x67, 0xe7, 0x17, 0x97, 0x57, 0xd7, 0x37, 0xb7, 0x77, 0xf7,
      0x0f, 0x8f, 0x4f, 0xcf, 0x2f, 0xaf, 0x6f, 0xef, 0x1f, 0x9f, 0x5f, 0xdf,
      0x3f, 0xbf, 0x7f, 0xff,
  };
  return table[x] >> (8 - n);
}

template <typename T> void dumph(T *data, uint32_t len, char sep = NULL) {
  for (int i = 0; i < len; i++)
    printf("%hhx%c", data[i], sep);
  std::cout << std::endl;
}
template <typename T>
void dumph(const char *msg, T *data, uint32_t len, char sep = NULL) {
  std::cout << msg;
  dumph<T>(data, len, sep);
}

class BitReader {

  inline static bool inited;

public:
  using bitbuf_type = uint8_t;
  inline static uint8_t bit_width = sizeof(bitbuf_type) * 8;
  inline static bitbuf_type bit_buffer;
  inline static uint8_t num_read;
  inline static std::vector<uint8_t> inbuf;
  inline static std::vector<uint8_t>::iterator in;

  static void init(std::vector<uint8_t> inV) {
    inbuf = inV;
    in = inbuf.begin();
    read_bitbuf();
    inited = true;
  }

  // returns false if we run out of data OR bad param
  // static bool read_bits(uint8_t n, uint8_t *res) {
  //  assert(inited);

  //  // No more bits left to read, cancel
  //  if (num_read == 0 && !read_bitbuf()) {
  //    return false;
  //  }

  //  *res = 0;
  //  uint8_t nbits = n;

  //  // If we carry over to then next number, need to remember how many digits
  //  // down we were before, so we dont overwrite digits we or'd in previously.
  //  // So we set this when we reset the bitbuffer to save that.
  //  uint8_t saved_shift = 0;

  //  // std::cout << "BBUF: " << std::bitset<8>(bit_buffer) << std::endl;
  //  std::cout << "REM: " << (int)num_read << std::endl;
  //  bool ret = true;
  //  if (n > sizeof(uint8_t) * 8) {
  //    std::cerr << "read_bits: n too big (" << (int)n << ")." << std::endl;
  //    ret = false;
  //    return ret;
  //  }

  //  uint8_t bits = 0;
  //  while (n-- > 0) {

  //    uint8_t shift = 8 - (nbits - n);
  //    // bits |= ((bit_buffer >> saved_shift) & (int(pow(2, shift))));
  //    bits |= ((bit_buffer >> saved_shift) & (int(pow(2, shift))));

  //    std::cout << std::bitset<8>(
  //                     ((bit_buffer >> saved_shift) & int(pow(2, shift))))
  //              << std::endl;
  //    // std::cout << std::bitset<8>(
  //    //                  ((bit_buffer >> saved_shift) & int(pow(2, shift))))
  //    //           << std::endl;
  //    num_read--;
  //    if (num_read != 0) {
  //      // bits <<= 1;
  //      continue;
  //    } else if (!read_bitbuf()) {
  //      ret = false;
  //      break;
  //    }

  //    // Assume we succesfully reset bitbuf, so can set the shift.

  //    saved_shift = 1;
  //    std::cout << "saved SHIFT: " << (int)saved_shift << std::endl;
  //  }

  //  bits >>= (bit_width - nbits);
  //  bit_buffer <<= 8 - num_read;
  //  //   if (bit_buffer == 0)
  //  //     read_bitbuf();
  //  std::cout << std::endl;
  //  std::cout << "BITS: " << std::bitset<8>(bits) << std::endl;
  //  std::cout << "BITS (as h): " << std::hex << (int)bits << std::dec
  //            << std::endl;
  //  *res = bits;
  //  return ret;
  //}

  static bool read_bits(uint8_t n, uint8_t *res) {
    assert(inited);

    uint8_t bits = 0;
    uint8_t nbits = n;
    uint8_t rem = 0;
    bool ret = true;

    // No more bits left to read, cancel
    if (num_read == 0 && !read_bitbuf()) {
      std::cout << "0: NREAD == 0" << std::endl;
      return false;
    }

    *res = 0;

    if (n > num_read)
      rem = n - num_read;

    uint8_t shift = 8 - n;
    uint8_t to_add = bit_buffer >> shift;
    bits |= to_add;

    num_read -= n + rem;
    nbits = n;

    if (rem) {
      if (!read_bitbuf()) {
        std::cout << "1: NREAD == 0" << std::endl;
        return false;
      }
      shift = 8 - rem;
      to_add = bit_buffer >> shift;
      bits |= to_add;

      num_read -= rem;
      nbits = rem;
    }

    bit_buffer <<= nbits;

    // std::cout << "Read BITS: " << std::bitset<8>(bits) << std::endl;

    *res = bits;

    return ret;
  }

  static bool read_bitbuf() {

    if (in == inbuf.end())
      return false;
    // std::cout << std::hex << (int)*in << std::dec << std::endl;
    bit_buffer = *in++;
    // Should be minimum bits or 8 bits if we still have bytes left to read
    num_read = 8;
    return true;
  }
};

// Class for writing bits
class BitWriter {

  inline static bool inited;

public:
  using bitbuf_type = uint8_t;
  inline static uint8_t bit_width = sizeof(bitbuf_type) * 8;
  inline static bitbuf_type bit_buffer;
  inline static uint8_t num_written;
  inline static std::vector<uint8_t> *inbuf;
  inline static std::vector<uint8_t>::iterator in;

  inline static std::vector<uint8_t> obuf;

  static void init(std::vector<uint8_t> *inV) {
    inbuf = inV;
    // Padding
    if (inbuf) {
      inbuf->push_back(0);
      in = inbuf->begin();
    }
    obuf.clear();
    bit_buffer = 0;
    num_written = 0;
    inited = true;
  }

  static void write_bits(uint8_t bits, uint8_t n) {
    assert(inited);

    if (n > bit_width) {
      std::cerr << "write_bits: n too big (" << (int)n << ")." << std::endl;
      return;
    }

    while (n-- > 0) {
      bit_buffer |= ((bits << (8 - (n + 1))) & 0x80) >> 7;
      // std::cout << std::bitset<8>(bit_buffer) << std::endl;
      num_written++;

      if (num_written == bit_width) {
        write_bitbuf();
      } else {

        bit_buffer <<= 1;
      }
    }
  }

  static void write_bitbuf() {
    uint8_t b0;
    b0 = bit_buffer;
    obuf.push_back(b0);
    num_written = 0;
    bit_buffer = 0;
  }

  static void write_in() {
    if (!inbuf)
      return;
    while (in != inbuf->end()) {
      std::cout << *in << std::endl;
      write_bits(*in, 8);
      in++;
    }
  }

  // flush out the bitbuffer, should write any unwritten bits
  static void flush_bits() {
    //  Just write whatever is left from the bitbuf
    //  Need to shift here as if we dont the remaining bits will be left on the
    //  low end of the bitbuffer, not ideal for reading.
    bit_buffer <<= (bit_width - num_written) - 1;
    write_bitbuf();
  }
};

// https://stackoverflow.com/questions/4954045/huffman-code-tables
//  n  bits   byte

// Basic huff table entry structure
struct huffc_tbl_entry {
  uint8_t n;
  unsigned char bits;
  unsigned char sym;
  // Does this reflect an actual entry or just a non-leaf - not written to file
  bool real;
};

// Wraps a huffs_tbl_entry, including a left and right branch used for encoding
struct huffc_encode_tbl_entry {
  bool isLeaf;
  struct huffc_tbl_entry entry;
  std::unique_ptr<struct huffc_encode_tbl_entry> left;
  std::unique_ptr<struct huffc_encode_tbl_entry> right;
  huffc_encode_tbl_entry() {
    memset(&entry, 0, sizeof(entry));
    // Prolly not necessary, i do it anyway
    left = 0;
    right = 0;
    isLeaf = true;
    entry.real = true;
  }
};

using huffc_tbl_type = std::vector<struct huffc_tbl_entry>;
using huffc_enc_tbl_type = std::unique_ptr<struct huffc_encode_tbl_entry>;

// TODO: Replace other uses of make_unique with dis
huffc_enc_tbl_type huffc_make_node() {
  return std::make_unique<struct huffc_encode_tbl_entry>(
      huffc_encode_tbl_entry());
}

using huffc_symtype = uint8_t;
using huffc_mtype = std::pair<huffc_symtype, uint32_t>;
using huffc_bits_type = uint8_t;

using huffc_bitvec_ele = std::pair<huffc_symtype, huffc_bits_type>;
using huffc_bitvec_type = std::vector<huffc_bitvec_ele>;

using huffc_sym_map = std::map<huffc_symtype, struct huffc_tbl_entry>;
using huffc_sym_index_map = std::unordered_map<huffc_symtype, uint32_t>;
// Compression stuff

// Increment an iterator safely
// Yes i know im lazy asf for just auto'ing everything
bool incIter(auto &iter, const auto &end) {
  if (iter++ == end)
    return false;
  return true;
}

#define INC_SYMC_ITER(x)                                                       \
  if (!incIter(huffc_symcount_iter, huffc_symcount_end))                       \
    return x;

// Make tree of 0 & 1s for the huffman entries
// Make sure we pass a reference since this func is designed to recurse, so we
// dont wanna copy the vector every time
// This func works but is rlly overkill, since we follow ever possible branch
// and combination. AKA, wayyy to many combos.
// bool huffc_add_bits_old(
//    huffc_bits_type curr_bits, huffc_bitvec_type &huffc_bits,
//    std::vector<huffc_mtype>::iterator huffc_symcount_iter,
//    const std::vector<huffc_mtype>::iterator &huffc_symcount_end) {
//
//  // First, add the current bits and symbol
//  huffc_symtype sym = huffc_symcount_iter->first;
//
//  if (!incIter(huffc_symcount_iter, huffc_symcount_end))
//    return false;
//
//  auto currp = std::pair<uint8_t, uint8_t>{sym, curr_bits};
//  huffc_bits.push_back(currp);
//
//  // Now, do same for left and right
//  // L
//  std::cout << "LHS" << std::endl;
//  huffc_add_bits_old(curr_bits << 1, huffc_bits, huffc_symcount_iter,
//                 huffc_symcount_end);
//
//  // R
//  std::cout << "RHS" << std::endl;
//  huffc_add_bits_old(curr_bits | 1, huffc_bits, huffc_symcount_iter,
//                 huffc_symcount_end);
//}

bool visualize_tree(struct huffc_encode_tbl_entry *root) {

  if (!root)
    return false;
  std::cout << "L: " << (char)(root->left ? root->left->entry.sym : 0)
            << "\tR: " << (char)(root->right ? root->right->entry.sym : 0)
            << std::endl;

  if (root->right)
    visualize_tree(root->right.get());
  if (root->left)
    visualize_tree(root->left.get());

  return false;
}

// For if we have a complete sequence of bits already
struct huffc_tbl_entry huffc_sym_lookup(huffc_enc_tbl_type &root,
                                        uint8_t bits) {

  struct huffc_encode_tbl_entry *curr = root.get();
  uint8_t bnum = 0, currb = 0;
  struct huffc_tbl_entry ent;

  std::cout << "Lookup BITS: " << std::bitset<8>(bits) << std::endl;

  while (curr) {

    ent = curr->entry;
    std::cout << "ENTRYBITS: " << std::bitset<8>(ent.bits)
              << "\t\tENTRYn: " << (int)ent.n << "\t\tENTRYSYM: " << ent.sym
              << std::endl;
    currb = (bits) & int(pow(2, bnum));

    // We found it
    // if (bnum == ent.n && bits == ent.bits)
    //  break;
    // if (bnum == ent.n && (bits & int(pow(2, bnum) - 1)) == ent.bits)

    if (currb) {

      // go right
      curr = curr->right.get();
    } else {
      curr = curr->left.get();
    }
    bnum++;
  }

  return ent;
}

// Returns entry on success, empty entry on failure
// If we need to lookup bits from the bit reader one at a time
struct huffc_tbl_entry huffc_sym_lookup(struct huffc_encode_tbl_entry *root) {

  auto curr = root;
  uint8_t bnum = 0, currb = 0;
  struct huffc_tbl_entry ent = {0};

  std::cout << "IN LOOKUP" << std::endl;

  while (curr) {

    ent = curr->entry;
    if (bnum == ent.n && ent.n != 0)
      break;
    // std::cout << "ENTRYBITS: " << std::bitset<8>(ent.bits)
    //           << "\t\tENTRYn: " << (int)ent.n << "\t\tENTRYSYM: " << ent.sym
    //           << std::endl;
    //   Reading bits one at a time because we dont yet have a way to "put back"
    //   bits we remove atm. Big TODO.
    if (!BitReader::read_bits(1, &currb)) {
      ent.n = 0;
      break;
    }
    // std::cout << (int)currb << " ";

    if (currb) {
      // go right
      curr = curr->right.get();
    } else {
      curr = curr->left.get();
    }
    bnum++;
  }

  return ent;
}

// This doesnt support codes longer than 8 bits, big todo.
// This is not done, and is fairly terrible, need to have proper ability to
// carry over bits too.
// uint8_t huffc_input = 0;
// uint8_t num_in_bits = 0;
// struct huffc_tbl_entry huffc_lookup(huffc_tbl_type tbl) {
//
//  if (num_in_bits == 0) {
//    BitReader::read_bits(8, &huffc_input);
//    num_in_bits = 8;
//  }
//
//  auto iter = tbl.begin();
//  while (iter != tbl.end()) {
//
//    huffc_tbl_entry e = *iter;
//
//    uint8_t realN = e.n;
//    uint8_t mask = ((int)pow(2, e.n) - 1) << (8 - (e.n - 1));
//
//    if ((huffc_input & mask) == e.bits) {
//      // Found, shift the input and other so we dont have to fuck around and
//      // break out.
//      huffc_input <<= e.n - 1;
//      num_in_bits -= e.n - 1;
//      break;
//    }
//
//    iter++;
//  }
//  return *iter;
//}

void huffc_add_element(struct huffc_tbl_entry ele,
                       struct huffc_encode_tbl_entry *root) {
  struct huffc_encode_tbl_entry *curr = root;
  uint32_t bnum = 0;
  std::cout << "BITTIES: " << std::bitset<8>(ele.bits) << std::endl;
  while (curr) {

    // Current bit
    uint8_t currb = ele.bits & int(pow(2, bnum));

    // We found our slot
    if (ele.n == bnum) {
      curr->entry = ele;
      return;
    }

    std::cout << "CURRB: " << std::bitset<8>(currb) << std::endl;
    if (currb) {
      std::cout << "R";
      // found branch, follow
      if (!curr->right) {
        curr->right = std::make_unique<struct huffc_encode_tbl_entry>(
            huffc_encode_tbl_entry());
      }
      curr = curr->right.get();
    } else {
      std::cout << "L";
      if (!curr->left) {
        // otherwise, make new branch
        curr->left = std::make_unique<struct huffc_encode_tbl_entry>(
            huffc_encode_tbl_entry());
      }
      curr = curr->left.get();
    }
    bnum++;
  }
  std::cout << std::endl;
}

// Returns trunk node of tree, or NULL if error
std::unique_ptr<struct huffc_encode_tbl_entry>
huffc_encode_tree(huffc_tbl_type &huffc_table) {

  auto root =
      std::make_unique<struct huffc_encode_tbl_entry>(huffc_encode_tbl_entry());
  auto it = huffc_table.begin();

  if (it == huffc_table.end()) {
    return NULL;
  }

  while (it != huffc_table.end()) {
    auto ele = *it++;
    huffc_add_element(ele, root.get());
  }

  return root;
}

bool is_pow2(uint32_t val) {
  // std::cout << val << " " << val - 1 << std::endl;
  return (val & (val - 1)) == 0;
}

// Big todo.
// Should recieve the iter for the map sorted in reverse order.
std::unique_ptr<struct huffc_encode_tbl_entry>
huffc_enc_tree(std::vector<huffc_mtype>::iterator huffc_symcount_iter,
               const std::vector<huffc_mtype>::iterator &huffc_symcount_end,
               huffc_tbl_type &tbl) {

  uint32_t curr_bits = 0;
  // Current num of bits
  uint32_t bit_ctr = 0;

  huffc_enc_tbl_type curr = NULL;

  // Stores the previous `curr`
  huffc_enc_tbl_type prev = NULL;

  while (huffc_symcount_iter != huffc_symcount_end) {

    struct huffc_tbl_entry e;

    huffc_mtype e0, e1;
    e1 = *huffc_symcount_iter;

    // Only move it if its not null, since idk what kinda fucky UB that would
    // cause
    if (curr != NULL) {
      // Move the ptr so prev refers to curr
      prev = std::move(curr);
    }

    if (prev == NULL) {
      INC_SYMC_ITER(curr)
      e0 = *huffc_symcount_iter;
    } else
      e0 = huffc_mtype({prev->entry.sym, prev->entry.n});

    // Nodes start out as non-leaf, so change fields to reflect that
    curr = huffc_make_node();
    curr->isLeaf = false;
    curr->entry.real = false;

    bool lt = (e0.second < e1.second);

    // Just straight up copy (move) over the entire entry
    if (lt && prev != NULL) {
      curr->left = std::move(prev);
      prev = NULL;
    } else {
      // Smaller goes to the left
      curr->left = huffc_make_node();
      curr->left->entry.bits = 0;
      curr->left->entry.sym = lt ? e0.first : e1.first;
      curr->left->entry.n = lt ? e0.second : e1.second;
    }

    // Store the entry in our table in increasing order of occurence.
    if (curr->left->entry.real)
      tbl.push_back(curr->left->entry);

    if (!lt && prev != NULL) {
      curr->right = std::move(prev);
      prev = NULL;
    } else {
      curr->right = huffc_make_node();
      curr->right->entry.bits = 1;
      curr->right->entry.sym = lt ? e1.first : e0.first;
      curr->right->entry.n = lt ? e1.second : e0.second;
    }

    if (curr->right->entry.real)
      tbl.push_back(curr->right->entry);

    std::cout << "N1: " << e0.first << "\tCOUNT: " << e0.second
              << "\t\t\tN2: " << e1.first << "\tCOUNT: " << e1.second
              << std::endl;
    // Now we need to count the total occurences of both symbols contained in
    // this node. Also need to carry this over so one of the next nodes we
    // process IS this node.
    curr->entry.n = e1.second + e0.second;

    // huffc_symcount_iter++;
    INC_SYMC_ITER(curr);
  }
  return curr;
}

void huffc_add_bits(
    huffc_tbl_type &huffc_table,
    std::vector<huffc_mtype>::iterator huffc_symcount_iter,
    const std::vector<huffc_mtype>::iterator &huffc_symcount_end) {

  uint32_t first_bits = 1;
  uint32_t curr_bits = 1;
  curr_bits <<= 1;
  // Count num of bits
  uint32_t bit_ctr = 1;
  // Only the first entry should have 1 in this slot
  huffc_bits_type sym_bits_curr = 1;

  while (huffc_symcount_iter != huffc_symcount_end) {

    huffc_symtype sym = huffc_symcount_iter->first;
    struct huffc_tbl_entry entry;
    entry.n = bit_ctr;
    entry.sym = sym;
    // entry.bits = reverse_byte(sym_bits_curr, entry.n);
    entry.bits = sym_bits_curr;

    std::cout << "BITS B4: " << std::bitset<8>(curr_bits) << std::endl;
    int ispow = (int)is_pow2(curr_bits);
    bit_ctr += ispow;
    // Keep the first entry clear and properly count the bits.
    sym_bits_curr = curr_bits++;
    if (sym_bits_curr & 1) {
      bit_ctr += (int)is_pow2(sym_bits_curr + 1);
      sym_bits_curr++;
      curr_bits++;
    }

    std::cout << "BITS AF: " << std::bitset<8>(curr_bits) << std::endl;

    huffc_table.push_back(entry);

    huffc_symcount_iter++;
  }
}

// Count symbol occurences, map to bits, creates table, tree and returns
huffc_enc_tbl_type huffc_count_syms(const huffc_symtype *data, uint32_t sz,
                                    huffc_tbl_type *tbl) {
  std::vector<huffc_mtype> huffc_symcount_flipped;

  { // anon namespace
    std::unordered_map<huffc_symtype, uint32_t> huffc_sym_count;

    // Count all occurences
    for (int i = 0; i < sz; i++) {
      huffc_sym_count[data[i]]++;
    }

    // Read map into vec, and sort that vec
    // For adding to the vector at back, basically just a push_back inserter
    auto inserter =
        std::back_inserter<std::vector<huffc_mtype>>(huffc_symcount_flipped);
    // Copy elements from map to vec
    std::copy(huffc_sym_count.begin(), huffc_sym_count.end(), inserter);
    // Sort the new vec
    std::sort(
        huffc_symcount_flipped.begin(), huffc_symcount_flipped.end(),
        [=](huffc_mtype a, huffc_mtype b) { return a.second < b.second; });
  } // anon namespace

  huffc_bits_type currBits = 0;
  // Map of sym -> bits
  // huffc_bitvec_type huffc_bits;
  std::vector<huffc_mtype>::iterator huffc_symcount_iter(
      huffc_symcount_flipped.begin());
  std::vector<huffc_mtype>::iterator huffc_symcount_end(
      huffc_symcount_flipped.end());

  huffc_tbl_type huffc_table;

  auto root =
      huffc_enc_tree(huffc_symcount_iter, huffc_symcount_end, huffc_table);

  std::cout << "HEAD: " << root->entry.sym << " N: " << root->entry.n
            << std::endl;
  visualize_tree(root.get());

  for (auto entry : huffc_table) {
    // TODO, figure this shit out
    std::cout << "char: " << entry.sym
              << " | bits: " << std::bitset<16>(entry.bits)
              << " | nbits: " << (int)entry.n << std::endl;
  }
  *tbl = huffc_table;
  return root;
}

int huffc_find_entry(huffc_symtype sym, huffc_tbl_type &tbl) {
  int idx = -1;
  for (int i = 0; i < tbl.size(); i++) {
    struct huffc_tbl_entry &ent = tbl[i];
    if (ent.sym != sym)
      continue;
    idx = i;
    break;
  }
  return idx;
}

std::vector<uint8_t> huffc_encode_data(const huffc_symtype *data, uint32_t sz,
                                       huffc_tbl_type &tbl) {

  BitWriter::init(0);

  // Now we need to write the bits into the buffer for each character
  huffc_sym_index_map sym_map;

  for (int i = 0; i < sz; i++) {
    // Havent seen, find index in vector
    huffc_symtype sym = data[i];
    if (!sym_map.contains(sym)) {
      int idx = huffc_find_entry(sym, tbl);
      assert(idx >= 0);
      sym_map[sym] = idx;
    } // Write bits

    auto &ent = tbl[sym_map[sym]];

    std::cout << "Symbol: " << sym
              << "\tWriting bits: " << std::bitset<8>(ent.bits)
              << "\tn: " << (int)ent.n << std::endl;

    // Write the bits reversed
    // BitWriter::write_bits(reverse_byte(ent.bits, ent.n), ent.n);
    BitWriter::write_bits(ent.bits, ent.n);
  }
  // Write any unwritten bits
  BitWriter::flush_bits();
  return BitWriter::obuf;
}

// Write the huffman table very very naively
std::vector<uint8_t> huffc_write_enc_and_tbl(huffc_tbl_type &tbl,
                                             std::vector<uint8_t> &enc) {

  std::vector<uint8_t> res;
  BitWriter::init(NULL);
  // First, we need to write a field denoting num of table entries.
  BitWriter::write_bits(tbl.size(), 8);
  for (auto entry : tbl) {
    BitWriter::write_bits(entry.n, 8);
    std::cout << "ENCODED N: " << (int)entry.n << std::endl;
    // BitWriter::write_bits(reverse_byte(entry.bits, entry.n), entry.n);
    BitWriter::write_bits(entry.bits, entry.n);
    std::cout << "ENCODED BITS: " << std::bitset<8>(entry.bits) << std::endl;
    BitWriter::write_bits(entry.sym, 8);
    std::cout << "ENCODED SYM: " << (char)entry.sym << std::endl;
  }

  for (auto byte : enc)
    BitWriter::write_bits(byte, 8);
  // Always remember to flush any unwritten bits
  BitWriter::flush_bits();
  res = BitWriter::obuf;
  return res;
}

// Write table and encoded data
std::vector<uint8_t> huffc_write_encoded(huffc_tbl_type &tbl,
                                         std::vector<uint8_t> &enc) {

  auto res = huffc_write_enc_and_tbl(tbl, enc);
  // res.insert(res.end(), enc.begin(), enc.end());

  return res;
}

huffc_tbl_type huffc_decode(std::vector<uint8_t> encoded) {

  huffc_tbl_type tbl;
  std::vector<uint8_t> decoded;
  BitReader::init(encoded);
  // First, get the table
  uint8_t num_entries = 0;
  BitReader::read_bits(8, &num_entries);
  std::cout << "Entries: " << (int)num_entries << std::endl;
  while (num_entries-- > 0) {
    huffc_tbl_entry e;
    BitReader::read_bits(8, &e.n);
    std::cout << "DECODED ENTRY N: " << (int)e.n << std::endl;
    BitReader::read_bits(e.n, &e.bits);
    std::cout << "DECODED ENTRY BITS: " << std::bitset<8>(e.bits) << std::endl;
    BitReader::read_bits(8, &e.sym);
    std::cout << "DECODED ENTRY SYM: " << (char)e.sym << std::endl;
    tbl.push_back(e);
  }

  // auto tree = huffc_encode_tree(tbl);

  // visualize_tree(tree.get());

  //// Now try to decode using this.
  // while (true) {
  //   auto ent = huffc_sym_lookup(tree.get());
  //   // No more entries OR failure in the lookup
  //   if (ent.n == 0) {
  //     std::cout << "N == 0" << std::endl;
  //     break;
  //   }
  //   std::cout << ent.sym << "\n\n";
  // }

  // Instead of that tree BS, just use the table we made

  uint8_t lookup_bits = 0, curr_bits = 0;
  uint8_t currN = 0;
  bool found = false;
  char curr_chr = 0;
  while (BitReader::read_bits(1, &lookup_bits)) {
    currN++;
    curr_bits <<= 1;
    curr_bits |= lookup_bits;

    // std::cout << "READ BITS: " << std::bitset<8>(curr_bits) << std::endl;

    for (int idx = tbl.size() - 1; idx >= 0; idx--) {
      auto e = tbl[idx];
      // std::cout << "E BITS: " << std::bitset<8>(e.bits) << std::endl;

      if (e.bits != curr_bits)
        continue;
      curr_chr = e.sym;
      found = true;
      break;
    }

    if (!found)
      continue;

    std::cout << curr_chr << std::endl;
    // If found, reset state, append to decoded and continue
    currN = 0;
    lookup_bits = 0;
    curr_chr = 0;
    curr_bits = 0;
    found = false;

    decoded.push_back(curr_chr);
  }

  exit(0);
  // std::cout << std::endl;
  return tbl;
}

uint8_t huffc_get_bit_enc(huffc_tbl_entry entry, huffc_enc_tbl_type &root,
                          uint8_t *nbits) {

  // Traverse the created tree, using the frequency to determine if we follow a
  // branch. We basically need to ensure, if our freq is less than both l and r
  // that we follow the NON-LEAF, because our symbol needs to be below that

  struct huffc_encode_tbl_entry *curr = root.get();

  std::cout << "Searching for: " << (char)entry.sym << std::endl;

  uint8_t bits = 0;
  // Count number of bits
  uint8_t ctr = 0;

  while (curr) {

    std::cout << ".";

    auto currEntry = curr->entry;

    bool lt = (entry.n < currEntry.n);
    if (!lt) {
      // We may have just found it
      std::cout << "Found(?). N: " << (int)curr->entry.n
                << "\tSYM: " << curr->entry.sym
                << "\tBITS: " << std::bitset<8>(bits) << std::endl;
      break;
    }

    // If not found, decide whether we go left or right.
    bool l;
    bool follow = !(curr->left->isLeaf) && curr->left->entry.n > entry.n;
    std::cout << "SYM: " << curr->entry.sym << " LEAF!" << std::endl;
    l = follow ||
        (curr->left->entry.n == entry.n && curr->left->entry.sym == entry.sym);
    bits <<= 1;
    if (l) {
      // No or with 0, cuz pointless
      curr = curr->left.get();
    } else {
      bits |= 1;
      curr = curr->right.get();
    }

    ctr++;
  }

  *nbits = ctr;
  return bits;
}

void huffc_get_bit_encodings(huffc_tbl_type &tbl, huffc_enc_tbl_type &root) {
  for (auto &entry : tbl) {
    // entry.n is meant to encode the number of bits we use, so set it as such.
    entry.bits = huffc_get_bit_enc(entry, root, &entry.n);
  }
}

int main() {

  huffc_symtype data[] =
      //    "AAAAAAACDGGHJHJJJ";
      //"BCAADDDCCACACAC";
      "BBBBCCCCCCCCaaaaaaaaaAAACVG";
  uint32_t dlen = sizeof(data) - 1;
  // First, need to count occurences of symbols, get their bit representation
  // according to this and store them as table entries.
  huffc_tbl_type huffc_tbl;
  huffc_enc_tbl_type root = huffc_count_syms(data, dlen, &huffc_tbl);

  std::vector<uint8_t> in(data, data + dlen);
  BitWriter::init(&in);
  BitWriter::write_bits(3, 3);
  BitWriter::write_in();
  BitWriter::flush_bits();

  std::cout << "dlen: " << dlen << std::endl;
  dumph(data, dlen, ' ');
  std::cout << "olen: " << BitWriter::obuf.size() << std::endl;
  dumph("\n", BitWriter::obuf.data(), BitWriter::obuf.size(), ' ');

  std::cout << std::endl;

  BitReader::init(BitWriter::obuf);
  uint8_t b = 0;
  BitReader::read_bits(3, &b);
  BitReader::read_bits(5, &b);
  BitReader::read_bits(3, &b);
  while (BitReader::read_bits(8, &b)) {
    std::cout << (char)b << " ";
  }
  std::cout << std::endl;

  huffc_get_bit_encodings(huffc_tbl, root);
  // Now that tree is made we need to encode our input.
  auto enc = huffc_encode_data(data, dlen, huffc_tbl);
  dumph("ENC: ", enc.data(), enc.size(), ' ');
  //    TODO now we need to encode the table and write it to the encoded data.
  //    Then commence the grand decoding.
  auto encoded_stream = huffc_write_encoded(huffc_tbl, enc);
  dumph("FULL STREAM: ", encoded_stream.data(), encoded_stream.size(), ' ');
  auto tbl = huffc_decode(encoded_stream);
}
