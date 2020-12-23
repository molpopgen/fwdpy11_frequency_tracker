#include <tuple>
#include <cstdint>
#include <map>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <fwdpy11/types/DiploidPopulation.hpp>
#include <fwdpy11/evolvets/SampleRecorder.hpp>

namespace py = pybind11;

struct FreqTracker
{
    // Origin time, position, effect size
    using key_type = std::tuple<std::uint32_t, double, double>;

    using trajectory = std::vector<std::uint32_t>;

    // Hashed associative container: key -> frequency trajectory
    // We are being a bit inefficient here, using a red-black
    // tree instead of a hashed container.  The reason is that C++
    // does not provide std::hash for tuples, even tuples of POD.
    // We could fix that, but won't here.
    using map_type = std::map<key_type, trajectory>;

    map_type trajectories;
    std::uint32_t burin_time;

    FreqTracker(std::uint32_t burnin) : trajectories{}, burin_time{burnin}
    {
    }

    inline void
    operator()(const fwdpy11::DiploidPopulation &pop,
               fwdpy11::SampleRecorder & /*unused*/)
    {
        if (pop.generation >= burin_time)
            {
                for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
                    {
                        if (pop.mcounts[i] > 0 && pop.mcounts[i] < 2 * pop.N)
                            {
                                auto key = std::make_tuple(pop.mutations[i].g,
                                                           pop.mutations[i].pos,
                                                           pop.mutations[i].s);
                                auto itr = trajectories.find(key);
                                if (itr == end(trajectories))
                                    {
                                        trajectories.emplace(std::move(key),
                                                             trajectory{pop.mcounts[i]});
                                    }
                                else
                                    {
                                        itr->second.push_back(pop.mcounts[i]);
                                    }
                            }
                    }
            }
    }
};

PYBIND11_MODULE(freqtracker, m)
{
    py::class_<FreqTracker>(m, "FreqTracker")
        .def(py::init<std::uint32_t>(), py::arg("burnin_time")) // This is __init__
        .def("__call__",
             &FreqTracker::operator()) // This makes it a valid Python callable
        .def_readonly(
            "trajectories",
            &FreqTracker::trajectories); // Get the data as a dict in Python.
}
