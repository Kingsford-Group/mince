#ifndef LIBRARY_FORMAT_HPP
#define LIBRARY_FORMAT_HPP

#include <ostream>

enum class ReadType { SINGLE_END, PAIRED_END };
enum class ReadOrientation { SAME, AWAY, TOWARD, NONE };
enum class ReadStrandedness { SA, AS, S, A, U };

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
    friend std::ostream& operator<<(std::ostream& os, LibraryFormat& lf);
};

#endif // LIBRARY_FORMAT_HPP
