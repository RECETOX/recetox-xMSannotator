// [[Rcpp::plugins("cpp11")]]

#include <vector>
#include <tuple>

#include <Rcpp.h>
using namespace Rcpp;


static bool almost_equal(double a, double b, double tolerance) noexcept {
    return std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * tolerance;
}

template <unsigned Index, typename Pack>
static NumericVector unpack(const std::vector<Pack>& packed) {
    return NumericVector::import_transform(std::begin(packed), std::end(packed), [] (const Pack& pack) {
        return std::get<Index>(pack);
    });
}

// [[Rcpp::export]]
DataFrame match_by_mass(DataFrame peaks, DataFrame adducts, DataFrame compounds, double tolerance) {
    const auto peak_id = as<NumericVector>(peaks["peak"]);
    const auto peak_mz = as<NumericVector>(peaks["mz"]);
    const auto adduct_id = as<CharacterVector>(adducts["adduct"]);
    const auto adduct_mz = as<NumericVector>(adducts["mass"]);
    const auto adduct_factor = as<NumericVector>(adducts["factor"]);
    const auto adduct_charge = as<NumericVector>(adducts["charge"]);
    const auto compound_id = as<NumericVector>(compounds["compound"]);
    const auto compound_mz = as<NumericVector>(compounds["monoisotopic_mass"]);

    auto matches = std::vector<std::tuple<R_xlen_t, R_xlen_t, R_xlen_t, double>>{};
    const auto n_peaks = peaks.nrows();
    const auto n_adducts = adducts.nrows();
    const auto n_compounds = compounds.nrows();

    for (R_xlen_t i = 0; i < n_adducts; ++i) {
        for (R_xlen_t j = 0; j < n_compounds; ++j) {
            const auto expected_mass = (adduct_factor[i] * compound_mz[j] + adduct_mz[i]) / adduct_charge[i];

            for (R_xlen_t k = 0; k < n_peaks; ++k) {
                if (almost_equal(peak_mz[k], expected_mass, tolerance)) {
                    matches.emplace_back(i, j, k, expected_mass);
                }
            }
        }
    }

    return DataFrame::create(
        Named("peak") = peak_id[unpack<2>(matches)],
        Named("adduct") = adduct_id[unpack<0>(matches)],
        Named("compound") = compound_id[unpack<1>(matches)],
        Named("expected_mass") = unpack<3>(matches)
    );
}
