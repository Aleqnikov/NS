#pragma once

#include <atomic>
#include <thread>
#include <vector>
#include <functional>

namespace parallel 
{

template <typename... Args>
void run(std::size_t numThreads, Args&&... args)
{
    std::atomic<int> curInd{0};
    std::vector<std::thread> threads;
    threads.reserve(numThreads);
    for (std::size_t t = 0; t < numThreads; ++t)
        threads.emplace_back(args..., std::ref(curInd), t);
    for (auto& t : threads)
        t.join();
}

}