#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// generic function for intersection
template <typename InputIterator1, typename InputIterator2>
inline double intersection(InputIterator1 begin1, InputIterator1 end1, 
                           InputIterator2 begin2) {
  
  // value to return
  double rval = 0;
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    
    // accumulate if appropirate
    if (d1 == 1 && d2 == 1)
      rval += 1;
  }
  return rval;  
}

// generic function for rowSum
template <typename InputIterator1, typename InputIterator2>
inline double allsum(InputIterator1 begin1,
                     InputIterator1 end1,
                     InputIterator2 begin2) {
  
  // value to return
  double rval = 0;
  
  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;
  
  // for each input item
  while (it1 != end1) {
    
    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;
    rval += d1 + d2;
  }
  return rval;  
}

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct JaccardDistance : public Worker {
  
  // input matrix to read from
  const RMatrix<double> mat;
  
  // output matrix to write to
  RMatrix<double> rmat;
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  JaccardDistance(const NumericMatrix mat, NumericMatrix rmat)
    : mat(mat), rmat(rmat) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (std::size_t j = 0; j < i; j++) {
        
        // rows we will operate on
        RMatrix<double>::Row row1 = mat.row(i);
        RMatrix<double>::Row row2 = mat.row(j);
        
        // calculate intersection and union
        double d1 = intersection(row1.begin(), row1.end(), row2.begin());
        double d2 = allsum(row1.begin(), row1.end(), row2.begin());
        
        // calculate Jaccard distance
        double d = d1/(d2 - d1);
        
        // calculate jaccard and write to output matrix
        rmat(i,j) = d;
        rmat(j,i) = d;
      }
    }
  }
};

//' Get Jaccard simmilarity for matrix rows
//'
//' @param mat Binary m x n matrix to be evaluated
//' @return Jaccard m x m simmilarity matrix
// [[Rcpp::export]]
NumericMatrix jaccard_fast_matrix(NumericMatrix mat) {
  
  // allocate the matrix we will return
  NumericMatrix rmat(mat.nrow(), mat.nrow());
  
  // create the worker
  JaccardDistance jaccardDistance(mat, rmat);
  
  // call it with parallelFor
  parallelFor(0, mat.nrow(), jaccardDistance);
  
  return rmat;
}
