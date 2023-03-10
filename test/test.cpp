#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "particle.hpp"

using Catch::Matchers::WithinRel;

TEST_CASE( "Particle sets mass correctly", "[particle]" ) {
    Particle p{3.14};
    REQUIRE_THAT( p.getMass(), WithinRel(3.14, 0.01) );
}
