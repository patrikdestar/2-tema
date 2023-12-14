#include <atomic>

class SpinBarrier
{
  private:
    std::atomic<int> count;             // Busy (not done) thread count
    std::atomic<int> cycle;             // Barrier iteration no.
    int maxc;                           // No. of participants (threads)
  public:
    SpinBarrier(int members)            /* Constructor */
        : count(members), cycle(0), maxc(members) { }

    void wait()                         /* Wait for others */
    {
        int c = cycle.load();
        if (--count == 0)               // If you are the last thread,
        {
            if (cycle.compare_exchange_weak(c,c+1))
                count = maxc;           // ...reset counter if needed...
            return;                     // ...and continue with code
        }
        while ((c == cycle) && (count != 0))
            std::this_thread::yield();  // Otherwise spin and wait
    }
};
