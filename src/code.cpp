#include <Rcpp.h>
using namespace Rcpp;

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
inline int rand_wrapper(const int n) { return floor(unif_rand()*n); }


//' Shuffles each column of a numeric matrix
//'
//' @param mat matrix to shuffle
//' @return matrix in which all the columns are shuffled (each column is shuffled separately)
// [[Rcpp::export]]
NumericMatrix shuffle_each_column(const NumericMatrix& mat){
	NumericMatrix m = Rcpp::clone(mat);    
    
    for (int i=0; i<m.ncol(); ++i){
        NumericMatrix::Column col = m( _ , i);
        
        std::random_shuffle(col.begin(), col.end(), rand_wrapper);        
    }      
    
    return(m);
}

// IntegerVector order_(const NumericVector& x) {
//   NumericVector sorted = clone(x).sort();
//   return match(sorted, x);
// }

// // [[Rcpp::export]]
// IntegerMatrix bounded_row_order(const NumericMatrix& mat, const int& k){    
//     int max_k = std::min(k, int(mat.ncol()));
    
//     IntegerMatrix m = IntegerMatrix(mat.nrow(), max_k);
    
//     for (int i=0; i<mat.nrow(); ++i){
//         NumericVector orig = mat(i, _);
//         IntegerMatrix::Row res = m(i, _);
        
//         IntegerVector ord = order_(orig);
                
//         std::copy( ord.begin(), ord.begin() + max_k, res.begin() ) ;

//         // Rcout << order_(orig) << "\n";
//     }      
    
//     return(m);
// }
