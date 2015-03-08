#ifndef __MINCE_OPTS_HPP__
#define __MINCE_OPTS_HPP__

#include "spdlog/spdlog.h"
#include <boost/filesystem.hpp>
#include <memory>

/**
 * This structure holds the common
 * options used by mince so that they
 * can easily be passed around together.
 */
struct MinceOpts {
    std::shared_ptr<spdlog::logger> jointLog{nullptr};
    std::shared_ptr<spdlog::logger> fileLog{nullptr};
    uint32_t numThreads;
};

#endif //__MINCE_OPTS_HPP__
