#include <TMB.hpp>
#include <algorithm>

// The code for splines is adopted from 'B-splines/QR decomposition'. 
// https://groups.google.com/g/tmb-users/c/AxtTSyONA7g/m/Z0GR0JdXAAAJ
template<class Type>
class SplineBasis {
public:
  int order,      /* order of the spline */
    ordm1,      /* order - 1 (3 for cubic splines) */
    nknots,     /* number of knots */
    curs,     /* current position in knots vector */
    boundary,   /* must have knots[(curs) <= x < knots(curs+1) */
    ncoef;      /* number of coefficients */
  /* except for the boundary case */
  vector<Type> ldel;    /* differences from knots on the left */
  vector<Type> rdel;  /* differences from knots on the right */
  vector<Type> knots; /* knot vector */
  vector<Type> coeff; /* coefficients */
  vector<Type> a;   /* scratch array */
  SplineBasis(int order = 4) : order(order) {
    ordm1 = order - 1;
    rdel = vector<Type>(ordm1);
    ldel = vector<Type>(ordm1);
    a = vector<Type>(order);
  }
  SplineBasis(vector<Type> knots, int order = 4) : order(order), knots(knots) {
    ordm1 = order - 1;
    nknots = knots.size();
    ncoef = nknots - order;
    rdel = vector<Type>(ordm1);
    ldel = vector<Type>(ordm1);
    a = vector<Type>(order);
  }
  int set_cursor(Type x)
  {
    int i;
    /* don't assume x's are sorted */
    curs = -1; /* Wall */
    boundary = 0;
    for (i = 0; i < nknots; i++) {
      if (knots(i) >= x) curs = i;
      if (knots(i) > x) break;
    }
    if (curs > nknots - order) {
      int lastLegit = nknots - order;
      if (x == knots(lastLegit)) {
        boundary = 1; curs = lastLegit;
      }
    }
    return curs;
  }
  void
  diff_table(Type x, int ndiff)
  {
    int i;
    for (i = 0; i < ndiff; i++) {
      rdel(i) = knots(curs + i) - x;
      ldel(i) = x - knots(curs - (i + 1));
    }
  }
  Type slow_evaluate(Type x, int nder)
  {
    int ti = curs, 
      lpt, apt, rpt, inner,
      outer = ordm1;
    if (boundary && nder == ordm1) { /* value is arbitrary */
      return Type(0);
    }
    while(nder--) {  // FIXME: divides by zero
      for(inner = outer, apt = 0, lpt = ti - outer; inner--; apt++, lpt++)
  a(apt) = Type(outer) * (a(apt + 1) - a(apt))/(knots(lpt + outer) - knots(lpt));
      outer--;
    }
    diff_table(x, outer);
    while(outer--)
      for(apt = 0, lpt = outer, rpt = 0, inner = outer + 1;
    inner--; lpt--, rpt++, apt++)
  // FIXME: divides by zero
  a(apt) = (a(apt + 1) * ldel(lpt) + a(apt) * rdel(rpt))/(rdel(rpt) + ldel(lpt));
    return a(0);
  }
  /* fast evaluation of basis functions */
  vector<Type> basis_funcs(Type x)
  {
    vector<Type> b(order);
    diff_table(x, ordm1);
    b(0) = Type(1);
    for (int j = 1; j <= ordm1; j++) {
      Type saved = Type(0);
      for (int r = 0; r < j; r++) { // do not divide by zero
  Type den = rdel(r) + ldel(j - 1 - r);
  if(den != Type(0)) {
    Type term = b(r)/den;
    b(r) = saved + rdel(r) * term;
    saved = ldel(j - 1 - r) * term;
  } else {
    if(r != Type(0) || rdel(r) != Type(0))
      b(r) = saved;
    saved = Type(0);
  }
      }
      b(j) = saved;
    }
    return b;
  }
  vector<Type> eval(Type x, int ders=0) {
    vector<Type> val(ncoef);
    val.setZero();
    set_cursor(x);
    int io = curs - order;
    if (io < 0 || io > nknots) {
      for (int j = 0; j < order; j++) {
  val(j+io) = Type(0); // R_NaN;
      }
    } else if (ders > 0) { /* slow method for derivatives */
      for(int i = 0; i < order; i++) {
  for(int j = 0; j < order; j++) a(j) = Type(0);
  a(i) = Type(1);
  val(i+io) = slow_evaluate(x, ders);
      }
    } else {    /* fast method for value */
      vector<Type> valtmp = basis_funcs(x);
      for (int i=0; i<valtmp.size(); i++)
        val(i+io)=valtmp(i);
    }
    return val;
  }
  matrix<Type> basis(vector<Type> x, int ders=0) {
    matrix<Type> mat(x.size(), ncoef);
    for (int i=0; i<x.size(); i++) {
      vector<Type> vec = eval(x(i), ders);
      for  (int j=0; j<vec.size(); j++)
  mat(i,j)=vec(j);
    }
    return mat;
  }
};
template<class Type>
class bs : public SplineBasis<Type> {
public:
  vector<Type> boundary_knots, interior_knots;
  int intercept, df;
  bs(vector<Type> boundary_knots, vector<Type> interior_knots, int intercept = 0) :
    SplineBasis<Type>(4), boundary_knots(boundary_knots), interior_knots(interior_knots),
    intercept(intercept) {
    df = intercept + 3 + interior_knots.size();
    this->nknots = interior_knots.size()+8;
    this->ncoef = this->nknots - this->order;
    this->knots = vector<Type>(this->nknots);
    for(int i=0; i<4;i++) {
      this->knots(i)=boundary_knots(0);
      this->knots(this->nknots-i-1)=boundary_knots(1);
    }
    if (interior_knots.size() > 0) 
      for(int i=0; i<interior_knots.size();i++) 
  this->knots(i+4)=interior_knots(i);
  }      
  vector<Type> eval(Type x, int ders=0) {
    vector<Type> vec;
    if (x<boundary_knots(0)) {
      Type k_pivot = Type(0.75)*boundary_knots(0)+Type(0.25)*interior_knots(0);
      Type delta = x - k_pivot;
      vec = bs<Type>::eval(k_pivot,0) +
        bs<Type>::eval(k_pivot,1)*delta +
        bs<Type>::eval(k_pivot,2)*delta*delta/Type(2) +
        bs<Type>::eval(k_pivot,3)*delta*delta*delta/Type(6);
    }
    else if (x>boundary_knots(1)) {
      Type k_pivot = Type(0.75)*boundary_knots(1)+Type(0.25)*interior_knots(interior_knots.size()-1);
      Type delta = x - k_pivot;
      vec = bs<Type>::eval(k_pivot,0) +
        bs<Type>::eval(k_pivot,1)*delta +
        bs<Type>::eval(k_pivot,2)*delta*delta/Type(2) +
        bs<Type>::eval(k_pivot,3)*delta*delta*delta/Type(6);
    }
    else  {
      vec = SplineBasis<Type>::eval(x, ders).segment(1-intercept,df);
    }
    return vec;
  }
  matrix<Type> basis(vector<Type> x, int ders=0) {
    matrix<Type> mat(x.size(), df);
    for (int i=0; i<x.size(); i++) {
      vector<Type> vec = bs<Type>::eval(x(i), ders);
      for  (int j=0; j<vec.size(); j++)
  mat(i,j)=vec(j);
    }
    return mat;
  }
};

template<class Type>
vector<Type> quantile(vector<Type> x, vector<Type> probs) {
  double alpha = 1.0, beta = 1.0;
  double n = x.size();
  std::vector<Type> x_copy = x;
  std::sort(x_copy.begin(), x_copy.end());
  vector<Type> inc_x = x_copy;
  vector<Type> res(probs.size());
  double fuzz = std::numeric_limits<double>::epsilon();
  for (int i = 0; i < probs.size(); i++) {
      // n * p + m
      double temp = asDouble(probs[i]);
      double nppm  = alpha + temp * (n + 1 - alpha - beta);
      double j = std::floor(nppm + fuzz);
      double h = nppm - j;
      long lj =  static_cast<long>(j);
      int pos = lj;
      Type res_i = (1 - h) * inc_x[pos - 1] + h * inc_x[pos];
      res[i] = res_i;
  }
  return res;
}

// help function: caclulate the Guass qaudarate integral for basis function
template<class Type> 
matrix<Type> guass_integral(vector<Type> quadwts_vec, matrix<Type> basis, int out_nquad) {
  matrix<Type> temp(quadwts_vec.size(), basis.cols());
  for (int i = 0; i < temp.cols(); i ++) {
    vector<Type> basiscoli = basis.col(i);
    temp.col(i) = quadwts_vec*basiscoli;
  }
  for (int i = 0; i < temp.cols(); i ++) {
    for (int j = 1; j < temp.rows(); j ++) {
      temp(j, i) = temp(j,i) + temp(j-1, i);
    }
  }
  matrix<Type> tempint(out_nquad, temp.cols());
  for (int i = 0; i < tempint.rows(); i ++) {
    vector<Type> tempinti = temp.row(i*2+1);
    tempint.row(i) = tempinti;
  }

  return tempint;
}

// help function: calculate the basis function according to the quantile
template<class Type>
matrix<Type> basis_fn(vector<Type> ixj, vector<Type> pos, vector<Type> quadwts_vec, int out_nquad) {
  Type min_ixj = min(ixj); Type max_ixj = max(ixj);
  Type ixj_center = 0.5 * (min_ixj + max_ixj); Type ixj_scale = (max_ixj - min_ixj);
  vector<Type> ixj_trans = (ixj - ixj_center)/ixj_scale;
  for (int i = 0; i < ixj_trans.size(); i ++) {
    ixj_trans[i] = 1/(1+exp(-ixj_trans[i]));
  }

  vector<Type> boundaryKnots(2);
  boundaryKnots[0] = min(ixj_trans); boundaryKnots[1] = max(ixj_trans);
  vector<Type> interiorKnots = quantile(ixj_trans, pos);

  bs<Type> basisj(boundaryKnots,interiorKnots);
  matrix<Type> basistemp = basisj.basis(ixj_trans);
  matrix<Type> basis_ixj = guass_integral(quadwts_vec, basistemp, out_nquad);
  for (int i = 0; i < basis_ixj.cols(); i++) {
    vector<Type> coli = basis_ixj.col(i);
    Type meani = coli.mean();
    basis_ixj.col(i) = coli - meani;
  }
  return basis_ixj;
}

// help function: calculate the inner integral in b prior
template<class Type>
matrix<Type> int_f(array<Type> a, matrix<Type> ix, vector<Type> pos, vector<Type> in_quadwts_vec, int out_nquad) {
  vector<int> a_dim = a.dim;
  matrix<Type> int_f(out_nquad, a_dim[1]);
  for (int i = 0;i < ix.cols(); i++) {
    matrix<Type> basistempi = basis_fn(vector<Type>(ix.col(i)), pos, in_quadwts_vec, out_nquad);
    matrix<Type> ai = a.col(i).matrix();
    int_f += basistempi*ai;
  }
  
  return int_f;
}


template <class Type>
Type objective_function<Type>::operator()()
{
  
  using namespace density;
  using namespace Eigen;
  
  
  // ------------------------------------------------------------------------ //
  // data
  // ------------------------------------------------------------------------ //
  DATA_MATRIX(y); 
  DATA_SCALAR(lambda);
  DATA_INTEGER(out_nquad);

  DATA_SPARSE_MATRIX(I0quadbasismat);
  DATA_SPARSE_MATRIX(I1quadbasismat);
  DATA_SPARSE_MATRIX(basismat);
  
  DATA_VECTOR(pos);
  DATA_VECTOR(eta);

  DATA_VECTOR(out_quadpts);
  DATA_VECTOR(out_quadwts);
  DATA_VECTOR(in_quadwts_vec);
  
  // ------------------------------------------------------------------------ //
  // parameters
  // ------------------------------------------------------------------------ //
  PARAMETER_VECTOR(logsigma);
  PARAMETER_MATRIX(b);
  PARAMETER_ARRAY(a);
  PARAMETER_MATRIX(theta);
  PARAMETER_VECTOR(a0);

  matrix<Type> xx = basismat * b;
  matrix<Type> ix = I1quadbasismat * b;
  matrix<Type> x = I0quadbasismat * b;

  
  int x_rows = x.rows();
  

  
  vector<Type> sigma = exp(logsigma);
  
  int m = a.rows();
  int nobs = y.rows();
  Type nll = 0.0;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from priors and data
  // ------------------------------------------------------------------------ //
  // likelihood
  for (int j = 0; j < y.cols(); j ++) {
    vector<Type> yj = y.col(j);
    vector<Type> xxj = xx.col(j);
    for(int i = 0; i < nobs; i++){
      nll -= dnorm(yj[i], xxj[i], sigma[j], true);
    }
  }

  // sigma prior
  for (int i = 0; i < logsigma.size(); i++) {
    nll -= logsigma[i];
  }

  // scale the matrix 
  for (int i = 0; i < x.cols(); i++) {
    vector<Type> coli = x.col(i);
    Type meani = coli.mean();
    x.col(i) = coli-meani;
  }

  // b prior  
  Type mean_out = out_quadpts.mean();
  matrix<Type> tempint = int_f(a, ix,pos, in_quadwts_vec, out_nquad);
  for (int j = 0; j < x.cols(); j ++) {
    vector<Type> xj = x.col(j);
    vector<Type> tempintj = tempint.col(j);
    vector<Type> a0j_int = a0[j]*(out_quadpts - mean_out);
    for (int i = 0; i < x_rows; i++) {
      nll += lambda/2*out_quadwts[i]*(xj[i] - a0j_int[i] - tempintj[i])*(xj[i] - a0j_int[i] - tempintj[i]);
    }
  }

  // theta prior
  for (int j = 0; j < theta.cols(); j ++) {
    vector<Type> thetaj = theta.col(j);
    for (int i = 0; i < theta.rows(); i++) {
      nll -= dbeta(thetaj[i], Type(1.0), Type(1.0), true);
    }
  }

  // a prior  
  for (int i = 0; i < x.cols(); i++) {
    for (int j = 0; j < x.cols(); j++) {
      vector<Type> aij = a.col(i).matrix().col(j);
      Type logpsij1 = m*log(eta[0]/sigma[i]) - (eta[0]/sigma[i])*sqrt((aij*aij).sum());
      Type logpsij2 = m*log(eta[1]/sigma[i]) - (eta[1]/sigma[i])*sqrt((aij*aij).sum());
      nll -= log((1-theta(j,i))*exp(logpsij1) + (theta(j,i))*exp(logpsij2));
    }
  }
  
  return nll;
}
