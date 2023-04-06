#pragma once

#include <stdint.h>

#if __linux__ != 0
#include <time.h>

static uint64_t timer_nsec() {
#if defined(CLOCK_MONOTONIC_RAW)
	const clockid_t clockid = CLOCK_MONOTONIC_RAW;

#else
	const clockid_t clockid = CLOCK_MONOTONIC;

#endif

	timespec t;
	clock_gettime(clockid, &t);

	return t.tv_sec * 1000000000UL + t.tv_nsec;
}

#elif _WIN64 != 0
#define NOMINMAX
#include <Windows.h>

static struct TimerBase {
	LARGE_INTEGER freq;
	TimerBase() {
		QueryPerformanceFrequency(&freq);
	}
} timerBase;

// the order of global initialisaitons is non-deterministic, do
// not use this routine in the ctors of globally-scoped objects
static uint64_t timer_nsec() {
	LARGE_INTEGER t;
	QueryPerformanceCounter(&t);

	return 1000000000ULL * t.QuadPart / timerBase.freq.QuadPart;
}

#elif __APPLE__ != 0
#include <mach/mach_time.h>

static struct TimerBase {
	mach_timebase_info_data_t tb;
	TimerBase() {
		mach_timebase_info(&tb);
	}
} timerBase;

// the order of global initialisaitons is non-deterministic, do
// not use this routine in the ctors of globally-scoped objects
static uint64_t timer_nsec() {
	const uint64_t t = mach_absolute_time();
	return t * timerBase.tb.numer / timerBase.tb.denom;
}

#endif

struct Timer {
	Timer()
		: start(timer_nsec()) {}

	template <typename T>
	static T toMs(T ns) {
		return T(ns / 1e6);
	}

	int64_t elapsedNs() const {
		const uint64_t now = timer_nsec();
		return now - start;
	}

	uint64_t start;
};

#include <cassert>

#define ASSERT_ALWAYS

#if defined(_MSC_VER)

#if defined(_DEBUG) || defined(ASSERT_ALWAYS)
#define MASSERT_ENABLED
#define mAssert(test) ( !!(test) ? (void)0 : __debugbreak() );
#else
#define mAssert(test) (void)0;
#endif
#else

#define mAssert(test) assert(test)

#endif

#include <vector>
#include <thread>
#include <mutex>

struct ThreadManager;

struct Task {
	/// Called before actually running the Task::run, always called on the same thread that is scheduling the task
	/// @param threadCount - number of threads executing this task
	virtual void onBeforeRun(int threadCount) {}

	/// The entry point for the thread executing this task
	/// @param threadIndex - 0 based index of the thread
	/// @param threadCount - number of threads executing this task
	virtual void run(int threadIndex, int threadCount) = 0;

	/// Run the task on a specific thread manager
	void runOn(ThreadManager &tm);

	virtual ~Task() {};
};

/// Non re-entrant task runner
struct ThreadManager {
	ThreadManager(const ThreadManager &) = delete;
	ThreadManager& operator=(const ThreadManager &) = delete;

	/// Initialize worker count without starting threads
	/// @param threadCount - the number of threads to spawn
	explicit ThreadManager(int threadCount)
		: count(threadCount)
	{}

	/// Start up all threads, must be called before @runThreads is called
	void start() {
		mAssert(count > 0 && "Task count must be positive");
		mAssert(threads.size() == 0 && "Already started");

		running = true;
		currentTask.resize(count, nullptr);

		threads.reserve(count);
		for (int c = 0; c < count; c++) {
			threads.emplace_back(&ThreadManager::threadBase, this, c);
		}
		// does not block to wait for threads to start, if they are delayed @runThreads will make sure to wait for them
	}

	/// Schedule the task to be run by the threads and return immediatelly
	/// This function could return before the threads have actually started running the task!
	void runThreadsNoWait(Task &task) {
		{
			std::lock_guard<std::mutex> lock(workMtx);
			mAssert(unlockedAllTasksDone() && "Already working");
			for (int c = 0; c < int(currentTask.size()); c++) {
				currentTask[c] = &task;
			}

			// inside the lock to ensure this is called before Task::run
			task.onBeforeRun(int(threads.size()));
		}

		workEvent.notify_all();
	}

	/// Start a task on all threads and wait for all of them to finish
	/// @param task - the task to run on all threads
	void runThreads(Task &task) {
		runThreadsNoWait(task);

		// block until task is complete
		if (!unlockedAllTasksDone()) {
			std::unique_lock<std::mutex> lock(workMtx);
			if (!unlockedAllTasksDone()) {
				workEvent.wait(lock, [this]() {
					return unlockedAllTasksDone();
				});
			}
		}
	}

	/// Blocking wait for all threads to exit, does not interrupt any running Task
	void stop() {
		mAssert(running && "Can't stop if not running");

		{
			std::lock_guard<std::mutex> lock(workMtx);
			running = false;
		}

		workEvent.notify_all();

		// joining threads will implicitly wait for all of them to finish current task
		for (int c = 0; c < int(threads.size()); c++) {
			threads[c].join();
		}

		threads.clear();
		currentTask.clear();
	}

	/// Get the number of worker threads
	int getThreadCount() const {
		return int(threads.size());
	}

private:

	/// The entry point for all of the threads
	/// @param threadIndex - the 0 based index of the thread
	void threadBase(volatile int threadIndex) {
		while (true) {

			Task *toExecute = nullptr;
			if (running) {
				std::unique_lock<std::mutex> lock(workMtx);

				if (running && currentTask[threadIndex] == nullptr) {
					workEvent.wait(lock, [this, threadIndex]() {
						return currentTask[threadIndex] != nullptr || !running;
					});
				}

				// just copy, and do not clear the value in @currentTask
				// it is used to signal the task is completed
				toExecute = currentTask[threadIndex];
			}

			if (!running) {
				return;
			}

			mAssert(toExecute);

			toExecute->run(threadIndex, int(threads.size()));
			{
				std::lock_guard<std::mutex> lock(workMtx);
				currentTask[threadIndex] = nullptr;
			}

			// Since start and finish share an event this must wake all threads
			// to make sure the caller is awoken and not only some other worker
			workEvent.notify_all();
		}
	}

	/// Check if all elements in @currentTask are nullptr
	/// @return true if at least one task is not nullptr, false otherwise
	bool unlockedAllTasksDone() const {
		for (int c = 0; c < currentTask.size(); c++) {
			if (currentTask[c]) {
				return false;
			}
		}
		return true;
	}

	int count = -1; ///< The number of threads
	std::vector<std::thread> threads; ///< The thread handles

	bool running = false; ///< Flag indicating if threads should quit

	/// The current task for each thread, always must be the same element
	/// Used to track if thread has finished working
	std::vector<Task *> currentTask;
	std::mutex workMtx; ///< Mutex protecting @currentTask and @running

	/// The event used to signal workers when new task is available
	/// Also used by workers to signal when task is finished
	std::condition_variable workEvent;
};


inline void Task::runOn(ThreadManager &tm) {
	tm.runThreads(*this);
}