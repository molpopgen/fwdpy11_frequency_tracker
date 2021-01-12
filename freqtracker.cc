#include <tuple>
#include <cstdint>
#include <unordered_map>
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

    // We record trajectories by nesting hashed associative containers.
    // The outer hash is on origin time.  The inner is on position.
    // Finally, the pair is effect size + the frequency trajectory.
    using map_type
        = std::unordered_map<std::uint32_t,
                             std::unordered_map<double, std::pair<double, trajectory>>>;

    map_type trajectories;
    std::uint32_t burnin_time;

    FreqTracker(std::uint32_t burnin) : trajectories{}, burnin_time{burnin}
    {
    }

    inline void
    operator()(const fwdpy11::DiploidPopulation &pop,
               fwdpy11::SampleRecorder & /*unused*/)
    {
        if (pop.generation >= burnin_time)
            {
                for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
                    {
                        if (pop.mcounts[i] > 0 && pop.mcounts[i] < 2 * pop.N)
                            {
                                auto origin = pop.mutations[i].g;
                                auto origin_itr = trajectories.find(origin);

                                if (origin_itr == end(trajectories))
                                    {
                                        trajectories[origin][pop.mutations[i].pos]
                                            = std::make_pair(pop.mutations[i].s,
                                                             trajectory{pop.mcounts[i]});
                                    }
                                else
                                    {
                                        auto pos_itr = origin_itr->second.find(
                                            pop.mutations[i].pos);
                                        if (pos_itr == end(origin_itr->second))
                                            {
                                                origin_itr->second[pop.mutations[i].pos]
                                                    = std::make_pair(
                                                        pop.mutations[i].s,
                                                        trajectory{pop.mcounts[i]});
                                            }
                                        else
                                            {
                                                pos_itr->second.second.push_back(
                                                    pop.mcounts[i]);
                                            }
                                    }
                            }
                    }
            }
    }

    // Make a copy of the data that Python will see
    // as a dict.
    // While we could rely on pybind11 to auto-convert
    // from a C++ map to a Python dict, doing so
    // would give 3 total copies in memory due to
    // temporaries.  So instead, we go straight
    // to a dict using pybind11's API.
    py::dict
    get_trajectories() const
    {
        py::dict rv;

        for (auto &&origin : trajectories)
            {
                for (auto &&position : origin.second)
                    {
                        // Make a Python tuple
                        auto key = py::make_tuple(origin.first, position.first,
                                                  position.second.first);
                        // Rely on auto-conversion of c++ vector
                        // to Python list, via pybind11/stl.h included above.
                        rv[key] = position.second.second;
                    }
            }

        return rv;
    }
};

PYBIND11_MODULE(freqtracker, m)
{
    py::class_<FreqTracker>(m, "FreqTracker")
        .def(py::init<std::uint32_t>(), py::arg("burnin_time")) // This is __init__
        .def("__call__",
             &FreqTracker::operator()) // This makes it a valid Python callable
        .def_property_readonly(
            "trajectories",
            &FreqTracker::get_trajectories); // Get the data as a dict in Python.
}
