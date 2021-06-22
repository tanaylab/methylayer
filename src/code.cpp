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

