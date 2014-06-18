/***************************************************************************
*                       Bit Stream File Class Implementation
*
*   File    : bitfile.cpp
*   Purpose : This file implements a simple class of I/O methods for files
*             that contain data in sizes that aren't integral bytes.  An
*             attempt was made to make the methods in this class analogous
*             to the methods provided to manipulate file streams.  The
*             methods contained in this class were created with compression
*             algorithms in mind, but may be suited to other applications.
*   Author  : Michael Dipperstein
*   Date    : July 20, 2004
*
****************************************************************************
*   UPDATES
*
*   $Id: bitfile.cpp,v 1.11 2009/07/23 03:57:56 michael Exp $
*   $Log: bitfile.cpp,v $
*   Revision 1.11  2009/07/23 03:57:56  michael
*   Zero out MSBs in value returned by GetBits when not retuning EOF.
*
*   Revision 1.10  2008/09/15 04:12:53  michael
*   Removed dead code.
*
*   Revision 1.9  2008/01/27 06:04:54  michael
*   Added  ByteAlign() and FlushOutput() methods.
*
*   Revision 1.8  2007/08/26 21:41:36  michael
*   All methods that don't modify the calling object have been made const to
*   increase functionality of const bit_array_c.
*
*   Changes required for LGPL v3.
*
*   Revision 1.7  2007/02/06 06:25:02  michael
*   Trim trailing spaces.
*
*   Revision 1.6  2006/06/03 18:57:59  michael
*   Corrected error discovered by anonymous in the destructor of writing
*   objects.  Underlying output stream was not being deleted.
*
*   Used spell checker to correct comments.
*
*   Revision 1.5  2006/02/10 04:30:47  michael
*   Applied fix for error discovered by Peter Husemann
*   <peter.husemann (at) cebitec (dot) uni-bielefeld (dot) de>.
*   When GetBit() reads a 0xFF byte, it would mistake it for EOF.
*
*   Revision 1.4  2005/12/10 05:20:01  michael
*   Added methods to get/put bits from/to integer types.
*
*   Revision 1.3  2005/06/23 04:39:06  michael
*   Convert from DOS end of line to Unix end of line
*
*   Revision 1.2  2005/06/23 04:33:07  michael
*   Prevent GetBits/PutBits from accessing an extra byte when given an
*   integral number of bytes.
*
*   Revision 1.1.1.1  2004/08/04 13:45:38  michael
*   bitfile class
*
*
****************************************************************************
*
* Bitfile: Bit Stream File I/O Class
* Copyright (C) 2004-2007 by Michael Dipperstein (mdipper@cs.ucsb.edu)
*
* This file is part of the bit file library.
*
* The bit file library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 3 of the
* License, or (at your option) any later version.
*
* The bit file library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
* General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
***************************************************************************/

/***************************************************************************
*                             INCLUDED FILES
***************************************************************************/
#include "bitfile.h"

using namespace std;

/***************************************************************************
*                            TYPE DEFINITIONS
***************************************************************************/

/* union used to test for endianess */
typedef union
{
    unsigned long word;
    unsigned char bytes[sizeof(unsigned long)];
} endian_test_t;

/***************************************************************************
*                                 METHODS
***************************************************************************/

/***************************************************************************
*   Method     : bit_file_c - default constructor
*   Description: This is the default bit_file_c constructor.  It
*                initializes stream pointers to NULL and clears the bit
*                buffer.
*   Parameters : None
*   Effects    : Initializes private members.
*   Returned   : None
***************************************************************************/
bit_file_c::bit_file_c(void)
{
    m_InStream = NULL;
    m_OutStream = NULL;
    m_BitBuffer = 0;
    m_BitCount = 0;
    m_Mode = BF_NO_MODE;

    /* test for endianess */
    endian_test_t endianTest;

    endianTest.word = 1;

    if (endianTest.bytes[0] == 1)
    {
        /* LSB is 1st byte (little endian)*/
        m_endian = BF_LITTLE_ENDIAN;
    }
    else if (endianTest.bytes[sizeof(unsigned long) - 1] == 1)
    {
        /* LSB is last byte (big endian)*/
        m_endian = BF_BIG_ENDIAN;
    }
    else
    {
        m_endian = BF_UNKNOWN_ENDIAN;
    }
}

/***************************************************************************
*   Method     : bit_file_c - constructor
*   Description: This is a bit_file_c constructor.  It opens an input or
*                output stream and clears the bit buffer.  An exception
*                will be thrown on error.
*   Parameters : fileName - NULL terminated string containing the name of
*                           the file to be opened.
*                mode - The mode of the file to be opened
*   Effects    : Initializes private members.  Creates and opens an input
*                or output stream.
*   Returned   : None
*   Exception  : "Error: Invalid File Type" - for unknown mode
*                "Error: Unable To Open File" - if stream cannot be opened
***************************************************************************/
bit_file_c::bit_file_c(const char *fileName, const BF_MODES mode)
{
    m_InStream = NULL;
    m_OutStream = NULL;
    m_BitBuffer = 0;
    m_BitCount = 0;

    switch (mode)
    {
        case BF_READ:
            m_InStream = new ifstream(fileName, ios::in | ios::binary);

            if (!m_InStream->good())
            {
                delete m_InStream;
                m_InStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }
            break;

        case BF_WRITE:
            m_OutStream = new ofstream(fileName, ios::out | ios::binary);

            if (!m_OutStream->good())
            {
                delete m_OutStream;
                m_OutStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }
            break;

        case BF_APPEND:
            m_OutStream =
                new ofstream(fileName, ios::out | ios::binary | ios::app);

            if (!m_OutStream->good())
            {
                delete m_OutStream;
                m_OutStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }
            break;

        default:
            throw("Error: Invalid File Type");
            break;
    }

    /* make sure we opened a file */
    if ((m_InStream == NULL) && (m_OutStream == NULL))
    {
        throw("Error: Unable To Open File");
    }

    /* test for endianess */
    endian_test_t endianTest;

    endianTest.word = 1;

    if (endianTest.bytes[0] == 1)
    {
        /* LSB is 1st byte (little endian)*/
        m_endian = BF_LITTLE_ENDIAN;
    }
    else if (endianTest.bytes[sizeof(unsigned long) - 1] == 1)
    {
        /* LSB is last byte (big endian)*/
        m_endian = BF_BIG_ENDIAN;
    }
    else
    {
        m_endian = BF_UNKNOWN_ENDIAN;
    }
}

/***************************************************************************
*   Method     : ~bit_file_c - destructor
*   Description: This is the bit_file_c destructor.  It closes and frees
*                any open file streams.  The bit buffer will be flushed
*                prior to closing an output stream.
*   Parameters : None
*   Effects    : Closes and frees open file streams.
*   Returned   : None
***************************************************************************/
bit_file_c::~bit_file_c(void)
{
    if (m_InStream != NULL)
    {
        m_InStream->close();
        delete m_InStream;
    }

    if (m_OutStream != NULL)
    {
        /* write out any unwritten bits */
        if (m_BitCount != 0)
        {
            m_BitBuffer <<= (8 - m_BitCount);
            m_OutStream->put(m_BitBuffer);
        }

        m_OutStream->close();
        delete m_OutStream;
    }
}

/***************************************************************************
*   Method     : Open
*   Description: This method opens an input or output stream and initializes
*                the bit buffer.  An exception will be thrown on error.
*   Parameters : fileName - NULL terminated string containing the name of
*                           the file to be opened.
*                mode - The mode of the file to be opened
*   Effects    : Creates and opens an input or output stream.
*   Returned   : None
*   Exception  : "Error: File Already Open" - if object has an open file
*                "Error: Invalid File Type" - for unknown mode
*                "Error: Unable To Open File" - if stream cannot be opened
***************************************************************************/
void bit_file_c::Open(const char *fileName, const BF_MODES mode)
{
    /* make sure file isn't already open */
    if ((m_InStream != NULL) || (m_OutStream != NULL))
    {
        throw("Error: File Already Open");
    }

    switch (mode)
    {
        case BF_READ:
            m_InStream = new ifstream(fileName, ios::in | ios::binary);

            if (!m_InStream->good())
            {
                delete m_InStream;
                m_InStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }

            m_BitBuffer = 0;
            m_BitCount = 0;
            break;

        case BF_WRITE:
            m_OutStream = new ofstream(fileName, ios::out | ios::binary);

            if (!m_OutStream->good())
            {
                delete m_OutStream;
                m_OutStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }

            m_BitBuffer = 0;
            m_BitCount = 0;
            break;

        case BF_APPEND:
            m_OutStream =
                new ofstream(fileName, ios::out | ios::binary | ios::app);

            if (!m_OutStream->good())
            {
                delete m_OutStream;
                m_OutStream = NULL;
            }
            else
            {
                m_Mode = mode;
            }

            m_BitBuffer = 0;
            m_BitCount = 0;
            break;

        default:
            throw("Error: Invalid File Type");
            break;
    }

    /* make sure we opened a file */
    if ((m_InStream == NULL) && (m_OutStream == NULL))
    {
        throw("Error: Unable To Open File");
    }
}

/***************************************************************************
*   Method     : Close
*   Description: This method closes and frees any open file streams.  The
*                bit buffer will be flushed prior to closing an output
*                stream.  All member variables are re-initialized.
*   Parameters : None
*   Effects    : Closes and frees open file streams.  Resets member
*                variables.
*   Returned   : None
***************************************************************************/
void bit_file_c::Close(void)
{
    if (m_InStream != NULL)
    {
        m_InStream->close();
        delete m_InStream;

        m_InStream = NULL;
        m_BitBuffer = 0;
        m_BitCount = 0;
        m_Mode = BF_NO_MODE;
    }

    if (m_OutStream != NULL)
    {
        /* write out any unwritten bits */
        if (m_BitCount != 0)
        {
            m_BitBuffer <<= (8 - m_BitCount);
            m_OutStream->put(m_BitBuffer);
        }

        m_OutStream->close();
        delete m_OutStream;

        m_OutStream = NULL;
        m_BitBuffer = 0;
        m_BitCount = 0;
        m_Mode = BF_NO_MODE;
    }
}

/***************************************************************************
*   Method     : ByteAlign
*   Description: This method aligns the bitfile to the nearest byte.  For
*                output files, this means writing out the bit buffer with
*                extra bits set to 0.  For input files, this means flushing
*                the bit buffer.
*   Parameters : None
*   Effects    : Flushes out the bit buffer.
*   Returned   : EOF if stream is NULL or write fails.  Writes return the
*                byte aligned contents of the bit buffer.  Reads returns
*                the unaligned contents of the bit buffer.
***************************************************************************/
int bit_file_c::ByteAlign(void)
{
    int returnValue;

    if ((BF_WRITE == m_Mode) || (BF_APPEND == m_Mode))
    {
        if (NULL == m_OutStream)
        {
            return(EOF);
        }
    }
    else
    {
        if (NULL == m_InStream)
        {
            return(EOF);
        }
    }

    returnValue = m_BitBuffer;

    if ((BF_WRITE == m_Mode) || (BF_APPEND == m_Mode))
    {
        /* write out any unwritten bits */
        if (m_BitCount != 0)
        {
            m_BitBuffer <<= 8 - (m_BitCount);
            m_OutStream->put(m_BitBuffer);  /* check for error */
        }
    }

    m_BitBuffer = 0;
    m_BitCount = 0;

    return (returnValue);
}

/***************************************************************************
*   Method     : FlushOutput
*   Description: This method flushes the output bit buffer.  This means
*                left justifying any pending bits, and filling spare bits
*                with the fill value.
*   Parameters : onesFill - non-zero if spare bits are filled with ones
*   Effects    : Flushes out the bit buffer, filling spare bits with ones
*                or zeros.
*   Returned   : EOF if stream is NULL or not writeable.  Otherwise, the
*                bit buffer value written. -1 if no data was written.
***************************************************************************/
int bit_file_c::FlushOutput(const unsigned char onesFill)
{
    int returnValue;

    if (NULL == m_OutStream)
    {
        return(EOF);
    }

    returnValue = -1;

    /* write out any unwritten bits */
    if (m_BitCount != 0)
    {
        m_BitBuffer <<= (8 - m_BitCount);

        if (onesFill)
        {
            m_BitBuffer |= (0xFF >> m_BitCount);
        }

        m_OutStream->put(m_BitBuffer);      /* check for error */
        returnValue = m_BitBuffer;
    }

    m_BitBuffer = 0;
    m_BitCount = 0;

    return (returnValue);
}

/***************************************************************************
*   Method     : GetChar
*   Description: This method returns the next byte from the input stream.
*   Parameters : None
*   Effects    : Reads next byte from file and updates buffer accordingly.
*   Returned   : EOF if a whole byte cannot be obtained.  Otherwise,
*                the character read.
***************************************************************************/
int bit_file_c::GetChar(void)
{
    int returnValue, tmp;

    if (m_InStream == NULL)
    {
        return EOF;
    }

    if (m_InStream->eof())
    {
        return EOF;
    }

    returnValue = m_InStream->get();

    if (m_BitCount == 0)
    {
        /* we can just get byte from file */
        return returnValue;
    }

    /* we have some buffered bits to return too */
    tmp = (returnValue >> m_BitCount);
    tmp |= m_BitBuffer << (8 - m_BitCount);

    /* put remaining in buffer. count shouldn't change. */
    m_BitBuffer = (char)returnValue;

    returnValue = tmp & 0xFF;

    return returnValue;
}

/***************************************************************************
*   Method     : PutChar
*   Description: This method writes the byte passed as a parameter to the
*                output stream.
*   Parameters : c - the character to be written
*   Effects    : Writes a byte to the file and updates buffer accordingly.
*   Returned   : On success, the character written, otherwise EOF.
***************************************************************************/
int bit_file_c::PutChar(const int c)
{
    int tmp;

    if (m_OutStream == NULL)
    {
        return EOF;
    }

    if (m_BitCount == 0)
    {
        /* we can just put byte from file */
        m_OutStream->put(c);
        return c;
    }

    /* figure out what to write */
    tmp = (c & 0xFF) >> m_BitCount;
    tmp = tmp | ((m_BitBuffer) << (8 - m_BitCount));

    m_OutStream->put((char)tmp);    /* check for error */

    /* put remaining in buffer. count shouldn't change. */
    m_BitBuffer = (char)c;

    return tmp;
}

/***************************************************************************
*   Method     : GetBit
*   Description: This method returns the next bit from the input stream.
*                The bit value returned is the msb in the bit buffer.
*   Parameters : None
*   Effects    : Reads next bit from bit buffer.  If the buffer is empty,
*                a new byte will be read from the input stream.
*   Returned   : 0 if bit == 0, 1 if bit == 1, and EOF if operation fails.
***************************************************************************/
int bit_file_c::GetBit(void)
{
    int returnValue;

    if (m_InStream == NULL)
    {
        return EOF;
    }

    if (m_BitCount == 0)
    {
        /* buffer is empty, read another character */
        if ((returnValue = m_InStream->get()) == EOF)
        {
            return EOF;         /* nothing left to read */
        }
        else
        {
            m_BitCount = 8;
            m_BitBuffer = returnValue;
        }
    }

    /* bit to return is msb in buffer */
    m_BitCount--;
    returnValue = m_BitBuffer >> m_BitCount;

    return (returnValue & 0x01);
}

/***************************************************************************
*   Method     : PutBit
*   Description: This method writes the bit passed as a parameter to the
*                output stream.
*   Parameters : c - the bit value to be written
*   Effects    : Writes a bit to the bit buffer.  If the buffer has a byte,
*                the buffer is written to the output stream and cleared.
*   Returned   : On success, the bit value written, otherwise EOF.
***************************************************************************/
int bit_file_c::PutBit(const int c)
{
    int returnValue = c;

    if (m_OutStream == NULL)
    {
        return EOF;
    }

    m_BitCount++;
    m_BitBuffer <<= 1;

    if (c != 0)
    {
        m_BitBuffer |= 1;
    }

    /* write bit buffer if we have 8 bits */
    if (m_BitCount == 8)
    {
        m_OutStream->put(m_BitBuffer);    /* check for error */

        /* reset buffer */
        m_BitCount = 0;
        m_BitBuffer = 0;
    }

    return returnValue;
}

/***************************************************************************
*   Method     : GetBits
*   Description: This method reads the specified number of bits from the
*                input stream and writes them to the requested memory
*                location (msb to lsb).
*   Parameters : bits - address to store bits read
*                count - number of bits to read
*   Effects    : Reads bits from the bit buffer and file stream.  The bit
*                buffer will be modified as necessary.
*   Returned   : EOF for failure, otherwise the number of bits read.  If
*                an EOF is reached before all the bits are read, bits
*                will contain every bit through the last complete byte.
***************************************************************************/
int bit_file_c::GetBits(void *bits, const unsigned int count)
{
    char *bytes, shifts;
    int offset, remaining, returnValue;

    if ((m_InStream == NULL) || (bits == NULL))
    {
        return EOF;
    }

    bytes = (char *)bits;

    offset = 0;
    remaining = count;

    /* read whole bytes */
    while (remaining >= 8)
    {
        returnValue = this->GetChar();

        if (returnValue == EOF)
        {
            return EOF;
        }

        bytes[offset] = (char)returnValue;
        remaining -= 8;
        offset++;
    }

    if (remaining != 0)
    {
        /* read remaining bits */
        shifts = 8 - remaining;

        while (remaining > 0)
        {
            returnValue = this->GetBit();

            if (returnValue == EOF)
            {
                return EOF;
            }

            bytes[offset] <<= 1;
            bytes[offset] |= (returnValue & 0x01);
            remaining--;
        }

        /* shift last bits into position */
        bytes[offset] <<= shifts;
    }

    return count;
}

/***************************************************************************
*   Method     : PutBits
*   Description: This method writes the specified number of bits from the
*                memory location passed as a parameter to the output
*                stream.   Bits are written msb to lsb.
*   Parameters : bits - pointer to bits to write
*                count - number of bits to write
*   Effects    : Writes bits to the bit buffer and file stream.  The bit
*                buffer will be modified as necessary.
*   Returned   : EOF for failure, otherwise the number of bits written.  If
*                an error occurs after a partial write, the partially
*                written bits will not be unwritten.
***************************************************************************/
int bit_file_c::PutBits(void *bits, const unsigned int count)
{
    char *bytes, tmp;
    int offset, remaining, returnValue;

    if ((m_OutStream == NULL) || (bits == NULL))
    {
        return EOF;
    }

    bytes = (char *)bits;

    offset = 0;
    remaining = count;

    /* write whole bytes */
    while (remaining >= 8)
    {
        returnValue = this->PutChar(bytes[offset]);

        if (returnValue == EOF)
        {
            return EOF;
        }

        remaining -= 8;
        offset++;
    }

    if (remaining != 0)
    {
        /* write remaining bits */
        tmp = bytes[offset];
        while (remaining > 0)
        {
            returnValue = this->PutBit(tmp & 0x80);

            if (returnValue == EOF)
            {
                return EOF;
            }

            tmp <<= 1;
            remaining--;
        }
    }

    return count;
}

/***************************************************************************
*   Method:    : GetBitsInt
*   Description: This method provides a machine independent layer that
*                allows a single call to stuff an arbitrary number of bits
*                read from the file stream into an integer type variable
*                (short, int, long, ...).
*   Parameters : bits - address to store bits read
*                count - number of bits to read
*                size - sizeof type containing "bits"
*   Effects    : Calls a method that reads bits from the bit buffer and
*                file stream.  The bit buffer will be modified as necessary.
*                the bits will be written to "bits" from least significant
*                byte to most significant byte.
*   Returned   : EOF for failure, otherwise the number of bits read by the
*                called method.  An error is thrown if the machine
*                endianess is unknown.
***************************************************************************/
int bit_file_c::GetBitsInt(void *bits, const unsigned int count,
    const size_t size)
{
    int returnValue;

    if ((m_InStream == NULL) || (bits == NULL))
    {
        return EOF;
    }

    if (m_endian == BF_LITTLE_ENDIAN)
    {
        returnValue = this->GetBitsLE(bits, count);
    }
    else if (m_endian == BF_BIG_ENDIAN)
    {
        returnValue = this->GetBitsBE(bits, count, size);
    }
    else
    {
        throw("Error: System Endianess Unknown");
    }

    return returnValue;
}

/***************************************************************************
*   Method     : GetBitsLE   (Little Endian)
*   Description: This method reads the specified number of bits from the
*                bit file and writes them to the requested memory location.
*                Bits are read LSB to MSB.
*   Parameters : bits - address to store bits read
*                count - number of bits to read
*   Effects    : Reads bits from the bit buffer and file stream.  The bit
*                buffer will be modified as necessary.  bits is treated as
*                a little endian integer of length >= (count/8) + 1.
*   Returned   : EOF for failure, otherwise the number of bits read.  If
*                an EOF is reached before all the bits are read, bits
*                will contain every bit through the last successful read.
***************************************************************************/
int bit_file_c::GetBitsLE(void *bits, const unsigned int count)
{
    char *bytes;
    int offset, remaining, returnValue;

    if ((m_InStream == NULL) || (bits == NULL))
    {
        return EOF;
    }

    bytes = (char *)bits;

    offset = 0;
    remaining = count;

    /* read whole bytes */
    while (remaining >= 8)
    {
        returnValue = this->GetChar();

        if (returnValue == EOF)
        {
            return EOF;
        }

        bytes[offset] = (char)returnValue;
        remaining -= 8;
        offset++;
    }

    if (remaining != 0)
    {
        /* read remaining bits */
        while (remaining > 0)
        {
            returnValue = this->GetBit();

            if (returnValue == EOF)
            {
                return EOF;
            }

            bytes[offset] <<= 1;
            bytes[offset] |= (returnValue & 0x01);
            remaining--;
        }
    }

    return count;
}

/***************************************************************************
*   Method     : GetBitsBE   (Big Endian)
*   Description: This method reads the specified number of bits from the
*                bit file and writes them to the requested memory location.
*                Bits are read LSB to MSB.
*   Parameters : bits - address to store bits read
*                count - number of bits to read
*                size - sizeof type containing "bits"
*   Effects    : Reads bits from the bit buffer and file stream.  The bit
*                buffer will be modified as necessary.  bits is treated as
*                a big endian integer of length size.
*   Returned   : EOF for failure, otherwise the number of bits read.  If
*                an EOF is reached before all the bits are read, bits
*                will contain every bit through the last successful read.
***************************************************************************/
int bit_file_c::GetBitsBE(void *bits, const unsigned int count,
    const size_t size)
{
    unsigned char *bytes;
    int offset, remaining, returnValue;

    if (count > (size * 8))
    {
        /* too many bits to read */
        return EOF;
    }

    bytes = (unsigned char *)bits;

    offset = size - 1;
    remaining = count;

    /* read whole bytes */
    while (remaining >= 8)
    {
        returnValue = this->GetChar();

        if (returnValue == EOF)
        {
            return EOF;
        }

        bytes[offset] = (unsigned char)returnValue;
        remaining -= 8;
        offset--;
    }

    if (remaining != 0)
    {
        /* read remaining bits */
        while (remaining > 0)
        {
            returnValue = this->GetBit();

            if (returnValue == EOF)
            {
                return EOF;
            }

            bytes[offset] <<= 1;
            bytes[offset] |= (returnValue & 0x01);
            remaining--;
        }
    }

    return count;
}

/***************************************************************************
*   Method     : PutBitsInt
*   Description: This method provides a machine independent layer that
*                allows a single function call to write an arbitrary number
*                of bits from an integer type variable (short, int, long,
*                ...) to the file stream.
*   Parameters : bits - pointer to bits to write
*                count - number of bits to write
*                size - sizeof type containing "bits"
*   Effects    : Calls a method that writes bits to the bit buffer and
*                file stream.  The bit buffer will be modified as necessary.
*                the bits will be written to the file stream from least
*                significant byte to most significant byte.
*   Returned   : EOF for failure, otherwise the number of bits written.  If
*                an error occurs after a partial write, the partially
*                written bits will not be unwritten.  An error is thrown if
*                the machine endianess is unknown.
***************************************************************************/
int bit_file_c::PutBitsInt(void *bits, const unsigned int count,
    const size_t size)
{
    int returnValue;

    if ((m_OutStream == NULL) || (bits == NULL))
    {
        return EOF;
    }

    if (m_endian == BF_LITTLE_ENDIAN)
    {
        returnValue = this->PutBitsLE(bits, count);
    }
    else if (m_endian == BF_BIG_ENDIAN)
    {
        returnValue = this->PutBitsBE(bits, count, size);
    }
    else
    {
        throw("Error: System Endianess Unknown");
    }

    return returnValue;
}

/***************************************************************************
*   Method     : PutBitsLE   (Little Endian)
*   Description: This method writes the specified number of bits from the
*                memory location passed as a parameter to the file stream
*                Bits are written LSB to MSB, assuming little endian order.
*   Parameters : bits - pointer to bits to write
*                count - number of bits to write
*   Effects    : Writes bits to the bit buffer and file stream.  The bit
*                buffer will be modified as necessary.  bits is treated as
*                a little endian integer of length >= (count/8) + 1.
*   Returned   : EOF for failure, otherwise the number of bits written.  If
*                an error occurs after a partial write, the partially
*                written bits will not be unwritten.
***************************************************************************/
int bit_file_c::PutBitsLE(void *bits, const unsigned int count)
{
    unsigned char *bytes, tmp;
    int offset, remaining, returnValue;

    bytes = (unsigned char *)bits;
    offset = 0;
    remaining = count;

    /* write whole bytes */
    while (remaining >= 8)
    {
        returnValue = this->PutChar(bytes[offset]);

        if (returnValue == EOF)
        {
            return EOF;
        }

        remaining -= 8;
        offset++;
    }

    if (remaining != 0)
    {
        /* write remaining bits */
        tmp = bytes[offset];
        tmp <<= (8 - remaining);

        while (remaining > 0)
        {
            returnValue = this->PutBit(tmp & 0x80);

            if (returnValue == EOF)
            {
                return EOF;
            }

            tmp <<= 1;
            remaining--;
        }
    }

    return count;
}

/***************************************************************************
*   Method     : PutBitsBE   (Big Endian)
*   Description: This method writes the specified number of bits from the
*                memory location passed as a parameter to the file stream
*                Bits are written LSB to MSB, assuming big endian order.
*   Parameters : bits - pointer to bits to write
*                count - number of bits to write
*   Effects    : Writes bits to the bit buffer and file stream.  The bit
*                buffer will be modified as necessary.  bits is treated as
*                a big endian integer of length size.
*   Returned   : EOF for failure, otherwise the number of bits written.  If
*                an error occurs after a partial write, the partially
*                written bits will not be unwritten.
***************************************************************************/
int bit_file_c::PutBitsBE(void *bits, const unsigned int count,
    const size_t size)
{
    unsigned char *bytes, tmp;
    int offset, remaining, returnValue;

    if (count > (size * 8))
    {
        /* too many bits to write */
        return EOF;
    }

    bytes = (unsigned char *)bits;
    offset = size - 1;
    remaining = count;

    /* write whole bytes */
    while (remaining >= 8)
    {
        returnValue = this->PutChar(bytes[offset]);

        if (returnValue == EOF)
        {
            return EOF;
        }

        remaining -= 8;
        offset--;
    }

    if (remaining != 0)
    {
        /* write remaining bits */
        tmp = bytes[offset];
        tmp <<= (8 - remaining);

        while (remaining > 0)
        {
            returnValue = this->PutBit(tmp & 0x80);

            if (returnValue == EOF)
            {
                return EOF;
            }

            tmp <<= 1;
            remaining--;
        }
    }

    return count;
}

/***************************************************************************
*   Method     : eof
*   Description: This method indicates whether or not the open file stream
*                is at the end of file.
*   Parameters : None
*   Effects    : None
*   Returned   : Returns true if the opened file stream is at an EOF.
*                Otherwise false is returned.
***************************************************************************/
bool bit_file_c::eof(void)
{
    if (m_InStream != NULL)
    {
        return (m_InStream->eof());
    }

    if (m_OutStream != NULL)
    {
        return (m_OutStream->eof());
    }

    /* return false for no file */
    return false;
}

/***************************************************************************
*   Method     : good
*   Description: This method is analogous to good for file streams.
*   Parameters : None
*   Effects    : None
*   Returned   : Returns good for the opened file stream.  False is
*                returned if there is no open file stream.
***************************************************************************/
bool bit_file_c::good(void)
{
    if (m_InStream != NULL)
    {
        return (m_InStream->good());
    }

    if (m_OutStream != NULL)
    {
        return (m_OutStream->good());
    }

    /* return false for no file */
    return false;
}

/***************************************************************************
*   Method     : bad
*   Description: This method is analogous to bad for file streams.
*   Parameters : None
*   Effects    : None
*   Returned   : Returns bad for the opened file stream.  False is
*                returned if there is no open file stream.
***************************************************************************/
bool bit_file_c::bad(void)
{
    if (m_InStream != NULL)
    {
        return (m_InStream->bad());
    }

    if (m_OutStream != NULL)
    {
        return (m_OutStream->bad());
    }

    /* return false for no file */
    return false;
}
