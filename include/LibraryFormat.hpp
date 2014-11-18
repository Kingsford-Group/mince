#ifndef LIBRARY_FORMAT_HPP
#define LIBRARY_FORMAT_HPP

#include <ostream>
#include <cstdint>

enum class ReadType : std::uint8_t { SINGLE_END = 0, PAIRED_END = 1 };
enum class ReadOrientation  : std::uint8_t { SAME = 0, AWAY = 1, TOWARD = 2, NONE = 3};
enum class ReadStrandedness  : std::uint8_t { SA = 0, AS = 1, S = 2, A = 3, U = 4 };

class LibraryFormat {
public:
    LibraryFormat(ReadType type_in, ReadOrientation orientation_in, ReadStrandedness strandedness_in);
    ReadType type;
    ReadOrientation orientation;
    ReadStrandedness strandedness;
    /**
     * Returns true if the specified library format is OK, false otherwise.
     */
    bool check();
    friend std::ostream& operator<<(std::ostream& os, const LibraryFormat& lf);
    // The maximum index that could be returned for a library type
    constexpr static uint8_t maxLibTypeID() {
        return 0 |
            ((0x01 & static_cast<uint8_t>(ReadType::PAIRED_END)) |
             (0xB & static_cast<uint8_t>(ReadOrientation::NONE)) << 1 |
             (0x6F & static_cast<uint8_t>(ReadStrandedness::U)) << 3);
    }

    inline static LibraryFormat formatFromID(uint8_t id) {
        ReadType rt;
        ReadOrientation ro;
        ReadStrandedness rs;
        switch (id & 0x01) {
            case 0:
                rt = ReadType::SINGLE_END;
                break;
            case 1:
                rt = ReadType::PAIRED_END;
                break;
        }
        switch ((id >> 1) & 0xB) {
            case static_cast<uint8_t>(ReadOrientation::SAME):
                ro = ReadOrientation::SAME;
                break;
            case static_cast<uint8_t>(ReadOrientation::AWAY):
                ro = ReadOrientation::AWAY;
                break;
            case static_cast<uint8_t>(ReadOrientation::TOWARD):
                ro = ReadOrientation::TOWARD;
                break;
            case static_cast<uint8_t>(ReadOrientation::NONE):
                ro = ReadOrientation::NONE;
                break;
        }

        switch ((id >> 3) & 0x6F) {
            case static_cast<uint8_t>(ReadStrandedness::SA):
                rs = ReadStrandedness::SA;
                break;
            case static_cast<uint8_t>(ReadStrandedness::AS):
                rs = ReadStrandedness::AS;
                break;
            case static_cast<uint8_t>(ReadStrandedness::S):
                rs = ReadStrandedness::S;
                break;
            case static_cast<uint8_t>(ReadStrandedness::A):
                rs = ReadStrandedness::A;
                break;
            case static_cast<uint8_t>(ReadStrandedness::U):
                rs = ReadStrandedness::U;
                break;
        }

        return LibraryFormat(rt, ro, rs);
    }

    // Assigns a unique ID to each potential library
    // type.  The IDs are such that 0 <= formatID(lib) < num possible formats
    inline uint8_t formatID() const {
        uint8_t id = 0;
        // mask: 00000001
        id |= (0x01 & static_cast<uint8_t>(type));
        // mask: 00000110
        id |= (0xB & static_cast<uint8_t>(orientation)) << 1;
        // mask: 00111000
        id |= (0x6F & static_cast<uint8_t>(strandedness)) << 3;
        return id;
    }
};

#endif // LIBRARY_FORMAT_HPP
