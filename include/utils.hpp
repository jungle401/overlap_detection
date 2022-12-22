#pragma once

#include <iostream>
#include <chrono>
#include <vector>

template <typename T>
void log(std::string interpretation, std::vector<T> v) {
  std::cout << interpretation << " :\t";
  for (auto& x : v) {
    std::cout << x << '\t';
  }
  std::cout << "\n";
}

template <typename T>
void log(T x, std::string interpretation) {
  std::cout << interpretation << " :\t";
  std::cout << x << '\n';
}

class Timer {
  public:
  std::chrono::duration<double> elapsed_seconds;
  std::chrono::time_point<std::chrono::steady_clock> time_point_start;
  std::chrono::time_point<std::chrono::steady_clock> time_point_end;
  std::string work;
  Timer() {}
  Timer(std::string work) {
		std::cout << " - [TIMER]: --> " << work << " ..." << std::endl;;
    this->time_point_start = std::chrono::steady_clock::now();
    this->work = work;
  }
  void start(std::string work) {
		std::cout << " - [TIMER]: --> " << work << " ..." << std::endl;;
    this->time_point_start = std::chrono::steady_clock::now();
    this->work = work;
  }
  void end() {
    this->time_point_end = std::chrono::steady_clock::now();
    this->elapsed_seconds = this->time_point_end - this->time_point_start;
    std::cout << " - [TIMER]: --| ";
    std::cout << this->work << ":\t";
    std::cout << elapsed_seconds.count() << ' ';
    std::cout << "seconds " << '\n';
  }
};
