/**
 * @file
 * @author  Giacomo Drago <giacomo@giacomodrago.com>
 * @version 1.0.1
 *
 *
 * @section LICENSE
 *
 * Copyright (c) 2013, Giacomo Drago <giacomo@giacomodrago.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *      This product includes fcmm, a software developed by Giacomo Drago.
 *      Website: http://projects.giacomodrago.com/fcmm
 * 4. Neither the name of Giacomo Drago nor the
 *    names of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GIACOMO DRAGO "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL GIACOMO DRAGO BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * @section DESCRIPTION
 *
 * This header file contains a template class implementing an almost lock-free concurrent memoization map.
 * Website: http://projects.giacomodrago.com/fcmm
 */

#ifndef FCMM_H_
#define FCMM_H_

// The noexcept specifier is unsupported in Visual Studio
#ifndef _MSC_VER
#define FCMM_NOEXCEPT noexcept
#else
#define FCMM_NOEXCEPT
#endif

#include <cstddef>
#include <cstdint>
#include <utility>
#include <algorithm>
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <type_traits>
#include <stdexcept>
#include <atomic>
#include <thread>

namespace fcmm {

/**
 * @brief Default maximum load factor
 */
static const float DEFAULT_MAX_LOAD_FACTOR = 0.75f;

/**
 * @brief Default maximum number of submaps
 */
static const std::size_t DEFAULT_MAX_NUM_SUBMAPS = 128;

/**
 * @brief The capacity of a new submap is calculated as the first prime number following the capacity of the
 * last submap multiplied by this constant
 */
static const std::size_t NEW_SUBMAPS_CAPACITY_MULTIPLIER = 8;

/**
 * @brief Minimum capacity of the first submap
 */
static const std::size_t FIRST_SUBMAP_MIN_CAPACITY = 65537;

/**
 * @brief The capacity of the first submap is calculated as:
 * <code>max(FIRST_SUBMAP_MIN_CAPACITY, nextPrime(<b>FIRST_SUBMAP_CAPACITY_MULTIPLIER</b> * estimatedNumEntries / maxLoadFactor))</code>
 */
static const float FIRST_SUBMAP_CAPACITY_MULTIPLIER = 1.03f;

/**
 * @brief Auxiliary function for nextPrime(). Returns `true` if `n` is prime, `false` otherwise.
 *
 * Adapted from http://stackoverflow.com/a/5694432/671092
 */
template<typename T>
static bool isPrime(const T& n, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr) {

    std::size_t divisor = 3;
    while (1) {
        const std::size_t quotient = n / divisor;
        if (quotient < divisor)
            return true;
        if (n == quotient * divisor)
            return false;
        divisor += 2;
    }

    return true;

}

/**
 * @brief Returns the smallest prime number greater than `n`
 *
 * Adapted from http://stackoverflow.com/a/5694432/671092
 */
template<typename T>
static T nextPrime(T n, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr) {

    if (n <= 2)
        return 2;

    if (n % 2 == 0)
        n++;

    while (!isPrime(n)) {
        n += 2;
    }

    return n;

}

/**
 * @brief This struct holds the statistics about a single submap of a @link Fcmm @endlink instance
 *
 * @see Stats
 */
struct SubmapStats {

    /**
     * @brief Capacity of the submap
     */
    std::size_t capacity;

    /**
     * @brief Number of valid buckets in the submap
     */
    std::size_t numValidBuckets;

    /**
     * @brief load factor of the submap (numValidBuckets / capacity)
     */
    float loadFactor;

};

/**
 * @brief This struct holds the statistics about a @link Fcmm @endlink instance
 *
 * @see Fcmm::getStats()
 */
struct Stats {

    /**
     * @brief Number of submaps in the map
     */
    std::size_t numSubmaps;

    /**
     * @brief Number of entries in the map
     */
    std::size_t numEntries;

    /**
     * @brief Statistics about each submap of the map
     * @see SubmapStats
     */
    std::vector<SubmapStats> submapsStats;

};

/**
 * @brief This hash function is the default for the `KeyHash2` template parameter of @link Fcmm @endlink.
 * It is only available for <a href="http://en.cppreference.com/w/cpp/types/is_integral">integral types</a>.
 */
template<typename Key, typename Enable = void>
class DefaultKeyHash2;

template<typename Key>
class DefaultKeyHash2<Key, typename std::enable_if<std::is_integral<Key>::value>::type> {
private:
    std::hash<Key> hash;
public:
    std::size_t operator()(const Key& key) const {
        return hash(~key);
    }
};

/**
 * @brief An almost lock-free concurrent hashmap, providing a limited set of functionalities.
 *
 *  - Only find() and insert() operations are supported: once an entry is inserted into the map, it cannot
 *    be erased nor updated. The filter() method provides a way to get a copy of the map containing only certain entries.
 *  - The map can be iterated over with an InputIterator (see begin() and end()) or a range-based for loop.
 *  - A deep copy of the map can be obtained with clone().
 *  - The presence of duplicate keys into the map is avoided but the total absence is not guaranteed.
 *    This is not an issue as long as it holds that: if (k<sub>1</sub>, v<sub>1</sub>) and (k<sub>2</sub>, v<sub>2</sub>)
 *    are entries and k<sub>1</sub> = k<sub>2</sub>, then v<sub>1</sub> = v<sub>2</sub>.
 *
 * Due to its features, this data structure is fit to be used for memoization in concurrent environments.
 *
 * @tparam  Key         the type of the key in each entry (default-constructible, copy-constructible or
 *                      move-constructible, copy-assignable or move-assignable)
 * @tparam  Value       the type of the value in each entry (default-constructible, copy-constructible or
 *                      move-constructible, copy-assignable or move-assignable)
 * @tparam  KeyHash1    the type of a function object that calculates the hash of the key;
 *                      it should have the same interface as
 *                      <a href="http://en.cppreference.com/w/cpp/utility/hash">`std::hash<T>`</a>
 * @tparam  KeyHash2    the type of another function object that calculates the hash of the key;
 *                      it should be <b>completely independent from `KeyHash1`</b>
 * @tparam  KeyEqual    the type of the function object that checks the equality of the two keys;
 *                      it should have the same interface as
 *                      <a href="http://en.cppreference.com/w/cpp/utility/functional/equal_to">`std::equal_to<T>`</a>
 */
template<
    typename Key,
    typename Value,
    typename KeyHash1 = std::hash<Key>,
    typename KeyHash2 = DefaultKeyHash2<Key>,
    typename KeyEqual = std::equal_to<Key>
>
class Fcmm {

    static_assert(std::is_default_constructible<Key>::value, "Key has to be default-constructible");
    static_assert(std::is_default_constructible<Value>::value, "Value has to be default-constructible");

    static_assert(std::is_copy_constructible<Key>::value || std::is_move_constructible<Key>::value,
                  "Key has to be copy-constructible or move-constructible, or both");
    static_assert(std::is_copy_constructible<Value>::value || std::is_move_constructible<Value>::value,
                  "Value has to be copy-constructible or move-constructible, or both");

    static_assert(std::is_copy_assignable<Key>::value || std::is_move_assignable<Key>::value,
                  "Key has to be copy-assignable or move-assignable, or both");
    static_assert(std::is_copy_assignable<Value>::value || std::is_move_assignable<Value>::value,
                  "Value has to be copy-assignable or move-assignable, or both");

public:

    /**
     * @brief The type of the key in each entry
     */
    typedef Key key_type;

    /**
     * @brief The type of the value in each entry
     */
    typedef Value mapped_type;

    /**
     * @brief An entry of the map
     */
    typedef std::pair<Key, Value> Entry;

    // Forward declaration
    class const_iterator;

private:

    /**
     * @brief First function object to calculate the hash code of a key
     */
    KeyHash1 keyHash1;

    /**
     * @brief Second function object to calculate the hash code of a key
     */
    KeyHash2 keyHash2;

    /**
     * @brief A bucket of the hashmap.
     *
     * A bucket can be in one of the following three states:
     *  - `EMPTY`: it does not contain an entry
     *  - `BUSY`: an entry is being written on it
     *  - `VALID`: it contains an entry
     */
    struct Bucket {

        enum class State { EMPTY, BUSY, VALID };

        std::atomic<State> state;
        Entry entry;

        Bucket() : state(State::EMPTY) {
        }

    };

    /**
     * @brief A submap is a collection of buckets
     */
    struct Submap {

        /**
         * @brief Function object to determine if two keys are equal
         */
        KeyEqual keyEqual;

        /**
         * @brief Buckets vector
         */
        std::vector<Bucket> buckets;

        /**
         * @brief Maximum load factor
         */
        float maxLoadFactor;

        /**
         * @brief Number of valid buckets
         */
        std::atomic<std::size_t> numValidBuckets;

        /**
         * @brief Constructor
         *
         * @param capacity       the capacity of this submap
         * @param maxLoadFactor  the maximum load factor of this submap
         */
        Submap(std::size_t capacity, float maxLoadFactor) :
                buckets(capacity),
                maxLoadFactor(maxLoadFactor),
                numValidBuckets(0) {
        }

        /**
         * @brief Returns the capacity of this submap
         */
        std::size_t getCapacity() const FCMM_NOEXCEPT {
            return buckets.size();
        }

        /**
         * @brief Returns a reference to a bucket of this submap
         *
         * @param index  the index of the bucket
         */
        Bucket& getBucket(std::size_t index) {
            return buckets[index];
        }

        /**
         * @brief Returns a const reference to a bucket of this submap
         *
         * @param index  the index of the bucket
         */
        const Bucket& getBucket(std::size_t index) const {
            return buckets[index];
        }

        /**
         * @brief Returns the number of entries in this submap
         */
        std::size_t getNumValidBuckets() const FCMM_NOEXCEPT {
            return numValidBuckets.load(std::memory_order_relaxed);
        }

        /**
         * @brief Increments the number of valid buckets by 1
         */
        void incrementNumValidBuckets() FCMM_NOEXCEPT {
            numValidBuckets.fetch_add(1, std::memory_order_relaxed);
        }

        /**
         * @brief Given the second hash of a key, calculates the corresponding
         * double hashing probe increment
         */
        std::size_t calculateProbeIncrement(std::size_t hash2) const FCMM_NOEXCEPT {
            const std::size_t modulus = getCapacity() - 1;
            return 1 + hash2 % modulus; // in [1, capacity - 1]
        }

        /**
         * @brief Searches for an entry having key equal to `key`
         *
         * @param key    the key of the entry to be found
         * @param hash1  the first hash of the key
         * @param hash2  the second hash of the key
         *
         * @return       a pair consisting of the index of the entry (if found)
         *               and a `bool` denoting whether the entry was found
         */
        std::pair<std::size_t, bool> find(const Key& key, std::size_t hash1, std::size_t hash2) const {

            const std::size_t startIndex = hash1 % getCapacity(); // initial position for probing
            const std::size_t probeIncrement = calculateProbeIncrement(hash2); // double hashing
            std::size_t index = startIndex; // current position for probing

            do {

                const Bucket& bucket = getBucket(index); // the current bucket being probed

                typename Bucket::State bucketState = bucket.state.load(std::memory_order_relaxed);

                if (bucketState == Bucket::State::VALID) {

                    std::atomic_thread_fence(std::memory_order_acquire); // memory fence

                    if (keyEqual(bucket.entry.first, key)) {
                        // the requested entry was found
                        return std::make_pair(index, true);
                    }

                } else if (bucketState == Bucket::State::EMPTY) {

                    // found a non-busy empty bucket: the requested entry is not present
                    return std::make_pair(0, false);

                }

                index = (index + probeIncrement) % getCapacity(); // move to the next bucket

            } while (index != startIndex);

            // scanned the whole submap: the requested entry is not present
            return std::make_pair(0, false);

        }

        /**
         * @brief Seeks the first valid bucket starting from `index` (inclusive)
         *
         * @param index  the index where the search should start;
         *               when the function returns, the variable contains the index where the valid bucket
         *               was found (if any)
         *
         * @return       `true` if a valid bucket was found, `false` if the end of the submap
         *               has been reached without finding a valid bucket
         */
        bool seek(std::size_t& index) const {

            while (index < getCapacity()) {

                const Bucket& bucket = getBucket(index); // the current bucket being probed

                if (bucket.state.load(std::memory_order_relaxed) == Bucket::State::VALID) {
                    std::atomic_thread_fence(std::memory_order_acquire); // memory fence
                    return true;
                }

                index++; // move to the next bucket

            }

            // reached the end of the buckets vector
            return false;

        }

        /**
         * @brief This exception is thrown by insert() if the new entry could not be inserted because the submap is full
         */
        struct FullSubmapException {
        };

        /**
         * @brief Inserts a new entry into the submap, if the submap doesn't already contain an entry with the same key.
         * The key will be moved, if possible.
         *
         * The value is calculated as needed by calling the `computeValue` function or functor.
         *
         * @param key                    the key of the entry to be inserted
         * @param hash1                  the first hash of the key
         * @param hash2                  the second hash of the key
         * @param computeValue           a function or functor that, given the key, calculates the corresponding value
         *
         * @return                       a pair consisting of the index of the entry (either inserted or preventing the insertion)
         *                               and a `bool` denoting whether the entry was inserted
         *
         * @tparam ComputeValueFunction  function or functor implementing `Value operator()(const Key&)`
         *
         * @throw FullSubmapException    thrown if the new entry could not be inserted because the submap is full
         *
         */
        template<typename KeyType, typename ComputeValueFunction>
        std::pair<std::size_t, bool> insert(KeyType&& key, std::size_t hash1, std::size_t hash2, ComputeValueFunction computeValue) {

            Value value = Value(); // to avoid "maybe-uninitialized" warnings
            bool valueComputed = false;

            const std::size_t startIndex = hash1 % getCapacity(); // initial position for probing
            const std::size_t probeIncrement = calculateProbeIncrement(hash2); // double hashing
            std::size_t index = startIndex; // current position for probing

            do {

                Bucket& bucket = getBucket(index); // the current bucket being probed

                typename Bucket::State bucketState = bucket.state.load(std::memory_order_relaxed);

                if (bucketState == Bucket::State::EMPTY) {

                    // since the bucket is (probably) empty, we will try to write the entry on it:
                    // let's compute the value of the entry (if it hasn't been computed yet)
                    if (!valueComputed) {
                        value = computeValue(key);
                        valueComputed = true;
                    }

                    // try to "lock" the bucket (without spinlocking)
                    if (bucket.state.compare_exchange_strong(bucketState, Bucket::State::BUSY, std::memory_order_relaxed)) {

                        // the bucket is now busy and this thread is the only one that can write on it

                        bucket.entry.first = std::move(key);
                        bucket.entry.second = std::move(value);
                        bucket.state.store(Bucket::State::VALID, std::memory_order_release); // mark the bucket as valid

                        incrementNumValidBuckets();

                        return std::make_pair(index, true);

                    }

                }

                // The following block cannot be turned into an else-if attatched to the previous if block, since the variable bucketState
                // may have been updated by compare_exchange_strong.
                // Moreover, if bucketState is different from VALID, we re-load a fresh value of the state variable and
                // check if it has become VALID in the meantime. This strategy reduces the presence of duplicates in the map.
                if (bucketState == Bucket::State::VALID || bucket.state.load(std::memory_order_relaxed) == Bucket::State::VALID) {

                    std::atomic_thread_fence(std::memory_order_acquire); // memory fence

                    if (keyEqual(bucket.entry.first, key)) { // does the key match?
                        // the key is already present in this submap: insertion failed
                        return std::make_pair(index, false);
                    }

                }

                index = (index + probeIncrement) % getCapacity(); // move to the next bucket

            } while (index != startIndex);

            // if the flow arrived here, then the submap is full: throw an exception that will be caught by Fcmm::insert()
            throw FullSubmapException();

        }

        /**
         * @brief Returns `true` if the submap is overloaded
         */
        bool isOverloaded() const FCMM_NOEXCEPT {
            return (float) getNumValidBuckets() / getCapacity() >= maxLoadFactor;
        }

        /**
         * @brief Returns statistics about this Fcmm::Submap instance.
         *
         * @see SubmapStats
         */
        SubmapStats getStats() const {

            SubmapStats stats;

            stats.capacity = getCapacity();
            stats.numValidBuckets = getNumValidBuckets();
            stats.loadFactor = (float) stats.numValidBuckets / stats.capacity;

            return stats;

        }

    };

    /**
     * @brief Maximum load factor
     */
    float maxLoadFactor;

    /**
     * @brief Number of submaps in this map
     */
    std::atomic<std::size_t> numSubmaps;

    /**
     * @brief Submaps pointers vector
     */
    std::vector<std::unique_ptr<Submap> > submaps;

    /**
     * @brief Number of entries in the map
     */
    std::atomic<std::size_t> numEntries;

    /**
     * @brief Atomic flag used by expand()
     */
    std::atomic_flag expanding;

    /**
     * @brief Returns the maximum number of submaps
     */
    std::size_t getMaxNumSubmaps() const FCMM_NOEXCEPT {
        return submaps.size();
    }

    /**
     * @brief Returns a reference to a submap pointer
     *
     * @param index  the index of the submap
     */
    std::unique_ptr<Submap>& getSubmap(std::size_t index) {
        return submaps[index];
    }

    /**
     * @brief Returns a const reference to a submap pointer
     *
     * @param index  the index of the submap
     */
    const std::unique_ptr<Submap>& getSubmap(std::size_t index) const {
        return submaps[index];
    }

    /**
     * @brief Returns the number of submaps
     */
    std::size_t getNumSubmaps() const FCMM_NOEXCEPT {
        return numSubmaps.load(std::memory_order_acquire);
    }

    /**
     * @brief Returns the index of the last submap
     */
    std::size_t getLastSubmapIndex() const FCMM_NOEXCEPT {
        return getNumSubmaps() - 1;
    }

    /**
     * @brief Increments the number of submaps by 1
     */
    void incrementNumSubmaps() FCMM_NOEXCEPT {
        numSubmaps.fetch_add(1, std::memory_order_release);
    }

    /**
     * @brief Increments the number of entries by 1
     */
    void incrementNumEntries() FCMM_NOEXCEPT {
        numEntries.fetch_add(1, std::memory_order_relaxed);
    }

    /**
     * @brief Expands the map by adding another submap
     *
     * @return  `true` if the map has been expanded, `false` otherwise
     */
    bool expand() {

        // spinlock on the expanding flag
        while (expanding.test_and_set(std::memory_order_acquire)) {
            std::this_thread::yield();
        }

        bool result = false;

        const std::size_t numSubmapsSnapshot = getNumSubmaps();

        // check if the maximum number of submaps has been reached
        if (numSubmapsSnapshot == getMaxNumSubmaps()) {
            throw std::runtime_error("Reached the maximum number of submaps: " + std::to_string(getMaxNumSubmaps()));
        }

        // get the last submap
        const std::size_t lastSubmapIndex = numSubmapsSnapshot - 1;
        const Submap& lastSubmap = *getSubmap(lastSubmapIndex);

        if (lastSubmap.isOverloaded()) { // re-check if the submap is overloaded
            // perform expansion
            const std::size_t newSubmapCapacity = nextPrime(lastSubmap.getCapacity() * NEW_SUBMAPS_CAPACITY_MULTIPLIER);
            getSubmap(lastSubmapIndex + 1).reset(new Submap(newSubmapCapacity, maxLoadFactor));
            incrementNumSubmaps();
            result = true;
        }

        // clear the expanding flag
        expanding.clear(std::memory_order_release);

        return result;

    }

    /**
     * @brief Searches for an entry having key equal to `key` across all the sumbaps in the range [0, lastSubmapIndex].
     *
     * @param key              the key of the entry to be found
     * @param hash1            the first hash of the key
     * @param hash2            the second hash of the key
     * @param lastSubmapIndex  the index of the last submap in which the entry has to be searched for
     *
     * @return                 a @link const_iterator @endlink to an entry having key `key`, or a past-the-end
     *                         const_iterator if no such entry is found
     */
    const_iterator findHelper(const Key& key, std::size_t hash1, std::size_t hash2, std::size_t lastSubmapIndex) const {

        for (long submapIndex = lastSubmapIndex; submapIndex >= 0; submapIndex--) { // scan each submap (from the last to the first)
            const Submap& submap = *getSubmap(submapIndex);
            const std::pair<std::size_t, bool> findResult = submap.find(key, hash1, hash2);
            if (findResult.second) { // the entry was found
                return const_iterator(this, submapIndex, findResult.first);
            }
        }

        return end(); // the entry was not found

    }

    /**
     * @brief Inserts a new entry into the map. The key will be moved, if possible.
     *
     * The value is calculated as needed by calling the `computeValue` function or functor.
     *
     * @param key                    the key of the entry to be inserted
     * @param hash1                  the first hash of the key
     * @param hash2                  the second hash of the key
     * @param computeValue           a function or functor that, given the key, calculates the corresponding value
     *
     * @tparam ComputeValueFunction  function or functor implementing `Value operator()(const Key&)`
     *
     * @return                       a pair consisting of a @link const_iterator @endlink to the inserted entry (or to the entry
     *                               that prevented the insertion) and a `bool` denoting whether the insertion took place
     */
    template<typename KeyType, typename ComputeValueFunction>
    std::pair<const_iterator, bool> insertHelper(KeyType&& key, std::size_t hash1, std::size_t hash2, ComputeValueFunction computeValue) {

        while (1) {

            const std::size_t lastSubmapIndex = getLastSubmapIndex();

            // check if the map (excl. the last submap) already contains a value for the key
            if (lastSubmapIndex > 0) {
                const const_iterator findIterator = findHelper(key, hash1, hash2, lastSubmapIndex - 1);
                if (findIterator != end()) { // the map already contains a value for the key
                    return std::make_pair(findIterator, false);
                }
            }

            Submap& lastSubmap = *getSubmap(lastSubmapIndex);

            if (lastSubmap.isOverloaded()) { // check if the submap is overloaded
                expand(); // expand the map
                continue; // restart the insertion process
            }

            try {
                const std::pair<std::size_t, bool> insertResult =
                        lastSubmap.insert(std::forward<KeyType>(key), hash1, hash2, computeValue);
                if (insertResult.second) {
                    incrementNumEntries();
                }
                const const_iterator insertIterator(this, lastSubmapIndex, insertResult.first);
                return std::make_pair(insertIterator, insertResult.second);
            } catch (typename Submap::FullSubmapException&) { // the submap is full
                expand(); // expand the map
                continue; // restart the insertion process
            }

        }

    }

public:

    /**
     * @brief Constructor
     *
     * @param estimatedNumEntries  an estimate of the number of entries this map will store
     * @param maxLoadFactor        the maximum load factor of each submap:
     *                             it should be a floating point number in the open interval (0, 1)
     * @param maxNumSubmaps        the maximum number of submaps that can be created (at least 1):
     *                             if this limit is exceeded, a `std::runtime_error` is thrown
     */
    Fcmm(std::size_t estimatedNumEntries = 0,
         float maxLoadFactor = DEFAULT_MAX_LOAD_FACTOR,
         std::size_t maxNumSubmaps = DEFAULT_MAX_NUM_SUBMAPS) :
            maxLoadFactor(maxLoadFactor),
            numSubmaps(1),
            submaps(maxNumSubmaps),
            numEntries(0) {

        // Not using ATOMIC_FLAG_INIT to workaround a Visual Studio bug
        expanding.clear();

        if (maxLoadFactor <= 0.0f || maxLoadFactor >= 1.0f) {
            throw std::logic_error("Invalid maximum load factor");
        }

        if (maxNumSubmaps < 1) {
            throw std::logic_error("Invalid maximum number of submaps");
        }

        // calculate the capacity of the first submap
        const std::size_t firstSubmapCapacity = std::max(
                    FIRST_SUBMAP_MIN_CAPACITY,
                    nextPrime((std::size_t) (FIRST_SUBMAP_CAPACITY_MULTIPLIER * estimatedNumEntries / maxLoadFactor)));

        // create the first submap
        getSubmap(0).reset(new Submap(firstSubmapCapacity, maxLoadFactor));

    }

    /**
     * @brief Searches for an entry having key equal to `key`
     *
     * @param key  the key of the entry to be found
     *
     * @return     an @link const_iterator @endlink to an entry having key `key`, or a past-the-end
     *             const_iterator if no such entry is found
     */
    const_iterator find(const Key& key) const {
        return findHelper(key, keyHash1(key), keyHash2(key), getLastSubmapIndex());
    }

    /**
     * @brief Returns a const reference to the value of an entry having key equal to `key`.
     * If no such entry exists, an exception of type `std::out_of_range` is thrown.
     *
     * @param key                the key of the entry to be found
     * @throw std::out_of_range  thrown if no entry exists having key equal to `key`
     */
    const Value& at(const Key& key) const {
        const const_iterator findIterator = find(key);
        if (findIterator == end()) {
            throw std::out_of_range("Entry not found");
        }
        return findIterator->second;
    }

    /**
     * @brief Returns a const reference to the value of an entry having key equal to `key`.
     * If no such entry exists, an exception of type `std::out_of_range` is thrown.
     *
     * @param key                the key of the entry to be found
     * @throw std::out_of_range  thrown if no entry exists having key equal to `key`
     */
    const Value& operator[](const Key& key) const {
        return at(key);
    }

    /**
     * @brief Inserts a new entry into the map.
     *
     * The value is calculated as needed by calling the `computeValue` function or functor.
     *
     * @param key                    the key of the entry to be inserted
     * @param computeValue           a function or functor that, given the key, calculates the corresponding value
     *
     * @tparam ComputeValueFunction  function or functor implementing `Value operator()(const Key&)`
     *
     * @return                       a pair consisting of a @link const_iterator @endlink to the inserted entry (or to the entry
     *                               that prevented the insertion) and a `bool` denoting whether the insertion took place
     */
    template<typename ComputeValueFunction>
    std::pair<const_iterator, bool> insert(const Key& key, ComputeValueFunction computeValue) {
        return insertHelper(key, keyHash1(key), keyHash2(key), computeValue);
    }

    /**
     * @brief Inserts a new entry into the map. The key will be moved, if possible.
     *
     * The value is calculated as needed by calling the `computeValue` function or functor.
     *
     * @param key                    the key of the entry to be inserted
     * @param computeValue           a function or functor that, given the key, calculates the corresponding value
     *
     * @tparam ComputeValueFunction  function or functor implementing `Value operator()(const Key&)`
     *
     * @return                       a pair consisting of a @link const_iterator @endlink to the inserted entry (or to the entry
     *                               that prevented the insertion) and a `bool` denoting whether the insertion took place
     */
    template<typename ComputeValueFunction>
    std::pair<const_iterator, bool> insert(Key&& key, ComputeValueFunction computeValue) {
        return insertHelper(key, keyHash1(key), keyHash2(key), computeValue);
    }

    /**
     * @brief Inserts a new entry into the map.
     *
     * @param entry  the entry to be inserted
     *
     * @return       a pair consisting of a @link const_iterator @endlink to the inserted entry (or to the entry
     *               that prevented the insertion) and a `bool` denoting whether the insertion took place
     *
     * @see insert(const Key&, ComputeValueFunction)
     */
    std::pair<const_iterator, bool> insert(const Entry& entry) {
        return insert(entry.first, [&entry](const Key&) -> const Value& { return entry.second; });
    }

    /**
     * @brief Inserts a new entry into the map. The members of the entry will be moved, if possible.
     *
     * @param entry  the entry to be inserted
     *
     * @return       a pair consisting of a @link const_iterator @endlink to the inserted entry (or to the entry
     *               that prevented the insertion) and a `bool` denoting whether the insertion took place
     *
     * @see insert(const Key&, ComputeValueFunction)
     */
    std::pair<const_iterator, bool> insert(Entry&& entry) {
        return insert(std::move(entry.first), [&entry](const Key&) -> Value&& { return std::move(entry.second); });
    }

    /**
     * @brief Inserts a new entry in the map. The new entry is constructed in place using `args` as the arguments
     * for the entry's constructor.
     *
     * @param args  arguments used to construct a new entry (key, value)
     *
     * @return      a pair consisting of a @link const_iterator @endlink to the inserted entry (or to the entry
     *              that prevented the insertion) and a `bool` denoting whether the insertion took place
     *
     * @see insert(const Key&, ComputeValueFunction)
     */
    template<typename... Args>
    std::pair<const_iterator, bool> emplace(Args&&... args) {
        return insert(Entry(std::forward<Args>(args)...));
    }

    /**
     * @brief Returns the number of entries in the map
     */
    std::size_t getNumEntries() const FCMM_NOEXCEPT {
        return numEntries.load(std::memory_order_relaxed);
    }

    /**
     * @brief Alias for getNumEntries()
     */
    std::size_t size() const FCMM_NOEXCEPT {
        return getNumEntries();
    }

    /**
     * @brief Returns `true` if the map has no elements, `false` otherwise
     */
    bool empty() const FCMM_NOEXCEPT {
        return getNumEntries() == 0;
    }

    /**
     * @brief Returns a @link const_iterator @endlink pointing to the first entry
     */
    const_iterator begin() const FCMM_NOEXCEPT {
        return const_iterator(this);
    }

    /**
     * @brief Returns a @link const_iterator @endlink pointing to the first entry
     */
    const_iterator cbegin() const FCMM_NOEXCEPT {
        return begin();
    }

    /**
     * @brief Returns a @link const_iterator @endlink pointing to the past-the-end entry
     */
    const_iterator end() const FCMM_NOEXCEPT {
        return const_iterator(this, true);
    }

    /**
     * @brief Returns a @link const_iterator @endlink pointing to the past-the-end entry
     */
    const_iterator cend() const FCMM_NOEXCEPT {
        return end();
    }

    /**
     * @brief Returns a pointer to a new map containing all the entries currently present in this map for
     * which `filterFunction(entry)` returns `true`.
     *
     * The new map is created via `new` and it is responsibility of the caller to `delete` it.
     *
     * @param filterFunction   a function or functor that, given an entry, returns `true` if
     *                         it should be kept, `false` if it should be filtered out
     *
     * @tparam FilterFunction  function or functor implementing `bool operator()(const Entry&)`
     *
     * @return                 the filtered map
     */
    template<typename FilterFunction>
    Fcmm* filter(FilterFunction filterFunction) const {

        Fcmm* map = new Fcmm(getNumEntries());

        for (const Entry& entry : *this) {
            if (filterFunction(entry)) {
                map->insert(entry);
            }
        }

        return map;

    }

    /**
     * @brief Returns a pointer to a new map containing all the entries currently present in this map.
     *
     * Actually, the number of entries in the cloned map can be lower than the number of entries in this map,
     * since duplicate entries (i.e. entries having the same key) are stripped out by the clone operation.
     *
     * The new map is created via `new` and it is responsibility of the caller to `delete` it.
     */

    Fcmm* clone() const {
        return filter([](const Entry&) { return true; });
    }

    /**
     * @brief Returns statistics about this @link Fcmm @endlink instance.
     *
     * The information returned by this method is useful for debugging and benchmarking.
     *
     * @see Stats
     */
    Stats getStats() const {

        Stats stats;

        stats.numEntries = getNumEntries();
        stats.numSubmaps = getNumSubmaps();
        for (std::size_t submapIndex = 0; submapIndex < stats.numSubmaps; submapIndex++) {
            const Submap& submap = *getSubmap(submapIndex);
            stats.submapsStats.push_back(submap.getStats());
        }

        return stats;

    }

    /**
     * @brief A const <a href="http://en.cppreference.com/w/cpp/concept/InputIterator">input iterator</a>
     * for iterating over the map. Iterators are never invalidated.
     */
    class const_iterator : public std::iterator<std::input_iterator_tag, const Entry> {

        friend class Fcmm;

    private:

        /**
         * @brief The @link Fcmm @endlink instance this iterator is iterating over
         */
        const Fcmm* map;

        /**
         * @brief The index of the submap containing the bucket currently pointed by this iterator
         */
        std::size_t submapIndex;

        /**
         * @brief The index of the bucket currently pointed by this iterator
         */
        std::size_t bucketIndex;

        /**
         * @brief This boolean is `true` if the iterator is pointing to the past-the-end entry of the map,
         * `false` otherwise
         */
        bool end;

        /**
         * @brief Returns the submap containing the bucket currently pointed by this iterator
         */
        const Submap& getSubmap() const {
            return *map->getSubmap(submapIndex);
        }

        /**
         * @brief Returns the bucket currently pointed by this iterator
         */
        const Bucket& getBucket() const {
            return getSubmap().getBucket(bucketIndex);
        }

        /**
         * @brief Returns the entry currently pointed by this iterator
         */
        const Entry& getEntry() const {
            return getBucket().entry;
        }

        /**
         * @brief Seeks the next valid bucket starting from `bucketIndex` (inclusive)
         */
        void seek() {

            while (!end) {

                if (getSubmap().seek(bucketIndex)) { // side-effect on bucketIndex
                    return;
                } else {
                    submapIndex++;
                    bucketIndex = 0;
                    if (submapIndex > map->getLastSubmapIndex()) {
                        end = true;
                        submapIndex = 0;
                    }
                }

            }

        }

        /**
         * @brief Increments `bucketIndex` and seeks the next valid bucket starting
         * from `bucketIndex` (inclusive)
         */
        void next() {
            bucketIndex++;
            seek();
        }

        /**
         * @brief Constructs a const iterator for a @link Fcmm @endlink instance
         *
         * @param map  the @link Fcmm @endlink instance the iterator shall iterate over
         * @param end  `true` if the iterator shall point to the past-the-end entry, `false` otherwise
         */
        const_iterator(const Fcmm* map, bool end = false) :
                map(map),
                submapIndex(0),
                bucketIndex(0),
                end(end) {
            seek();
        }

        /**
         * @brief Constructs a const iterator for a @link Fcmm @endlink instance
         *
         * @param map          the @link Fcmm @endlink instance the iterator shall iterate over
         * @param submapIndex  the index of the submap containing the bucket the iterator shall point to
         * @param bucketIndex  the index of the bucket the iterator shall point to
         */
        const_iterator(const Fcmm* map, std::size_t submapIndex, std::size_t bucketIndex) :
                map(map),
                submapIndex(submapIndex),
                bucketIndex(bucketIndex),
                end(false) {
        }

    public:

        /**
         * @brief Equality operator
         */
        bool operator==(const const_iterator& other) const {
            return map == other.map &&
                    ((end && other.end) || (!end && !other.end && submapIndex == other.submapIndex && bucketIndex == other.bucketIndex));
        }

        /**
         * @brief Inequality operator
         */
        bool operator!=(const const_iterator& other) const {
            return !(*this == other);
        }

        /**
         * @brief Pre-increment operator
         */
        const_iterator& operator++(void) {
            next();
            return *this;
        }

        /**
         * @brief Post-increment operator
         */
        const_iterator operator++(/* dummy */ int) {
            const_iterator old(*this);
            ++(*this);
            return old;
        }

        /**
         * @brief Dereference operator
         */
        const Entry& operator*() const {
            return getEntry();
        }

        /**
         * @brief Arrow operator
         */
        const Entry* operator->() const {
            return &getEntry();
        }

        /**
         * @brief Swap function
         */
        friend void swap(const_iterator& it1, const_iterator& it2) {
            std::swap(it1.map, it2.map);
            std::swap(it1.submapIndex, it2.submapIndex);
            std::swap(it1.bucketIndex, it2.bucketIndex);
            std::swap(it1.end, it2.end);
        }

    };

    /**
     * @brief The object is not copy-constructible: use clone() instead
     */
    Fcmm(const Fcmm&) = delete;

    /**
     * @brief The object cannot be moved
     */
    Fcmm(Fcmm&&) = delete;

};

} // namespace fcmm

#undef FCMM_NOEXCEPT

#endif // FCMM_H_
