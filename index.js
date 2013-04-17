
/* ~ App deps ~ */

var _ = require('underscore');

module.exports = Polyfit;

/**
 * Polyfit constructor
 */

function Polyfit(x, y) {

  // Make sure we return an instance
  if (!(this instanceof Polyfit)) return new Polyfit(x, y);

  // Array validity check
  if (!_.isArray(x) || !_.isArray(y)) {
    throw new Error('x and y must be arrays');
  }

  // Make sure we have equal lengths
  if (x.length != y.length) {
    throw new Error('x and y must have the same length');
  }

  this.matrix = [];

  // Zip the x and y values up into a matrix
  for (var i = 0, len = x.length; i < len; i++) {
    this.matrix[ i ] = { x: x[ i ], y: y[ i ] };
  }

}

/**
 * Perform gauss-jordan division
 */

Polyfit.prototype.gaussJordanDivide = function(matrix, row, col, numCols) {

  for (var i = col + 1; i < numCols; i++) {
    matrix[ row ][ i ] /= matrix[ row ][ col ];
  }

  matrix[ row ][ col ] = 1;

  return matrix;

};

/**
 * Perform gauss-jordan elimination
 */

Polyfit.prototype.gaussJordanEliminate = function(matrix, row, column, numRows, numCols) {

  for (var i = 0; i < numRows; i++) {
    if (i != row && matrix[ i ][ column ] != 0) {
      for (var j = column + 1; j < numCols; j++) {
        matrix[ i ][ j ] -= matrix[ i ][ column ] * matrix[ row ][ j ];
      }
      matrix[ i ][ column ] = 0;
    }
  }

  return this.matrix;

};

/**
 * Perform gauss-jordan echelon method
 */

Polyfit.prototype.gaussJordanEchelonize = function(matrix) {

  var rows = matrix.length
    , cols = matrix[ 0 ].length
    , i = 0
    , j = 0
    , k
    , swap;

  while (i < rows && j < cols) {
    k = i;
    // Look for non-zero entries in col j at or below row i
    while (k < rows && matrix[ k ][ j ] == 0) {
      k++;
    }
    // If an entry is found at row k
    if (k < rows) {
      // If k is not i, then swap row i with row k
      if (k != i) {
        swap = matrix[ i ];
        matrix[ i ] = matrix[ k ];
        matrix[ k ] = swap;
      }
      // If matrix[i][j] is != 1, divide row i by matrix[i][j]
      if (matrix[ i ][ j ] != 1) {
        this.gaussJordanDivide(matrix, i, j, cols);
      }
      // Eliminate all other non-zero entries
      this.gaussJordanEliminate(matrix, i, j, rows, cols);
      i++;
    }
    j++;
  }

  return matrix;

};

/**
 * Perform regression
 */

Polyfit.prototype.regress = function(x, terms) {

  var a = 0
    , exp = 0;

  for (var i = 0, len = terms.length; i < len; i++) {
    a += terms[ i ] * Math.pow(x, exp++);
  }

  return a;

};

/**
 * Compute correlation coefficient
 */

Polyfit.prototype.correlationCoefficient = function(terms) {

  var r = 0
    , n = this.matrix.length
    , sx = 0
    , sx2 = 0
    , sy = 0
    , sy2 = 0
    , sxy = 0
    , x
    , y;

  for (var i = 0; i < n; i++) {
    var pr = this.matrix[ i ];
    x = this.regress(pr.x, terms);
    y = pr.y;
    sx += x;
    sy += y;
    sxy += x * y;
    sx2 += x * x;
    sy2 += y * y;
  }

  var div = Math.sqrt((sx2 - (sx * sx) / n) * (sy2 - (sy * sy) / n));

  if (div != 0) {
    r = Math.pow((sxy - (sx * sy) / n) / div, 2);
  }

  return r;

};

/**
 * Run standard error function
 */

Polyfit.prototype.standardError = function(terms) {

  var r = 0
    , n = this.matrix.length;

  if (n > 2) {
    var a = 0;
    for (var i = 0; i < n; i++) {
      var pr = this.matrix[ i ];
      a += Math.pow((this.regress(pr.x, terms) - pr.y), 2);
    }
    r = Math.sqrt(a / (n - 2));
  }

  return r;

};

/**
 * Compute coefficients for given data matrix
 */

Polyfit.prototype.computeCoefficients = function(p) {

  var n = this.matrix.length
    , r
    , c
    , rs = 2 * (++p) - 1;

  var m = [];

  // Initialize array with 0 values
  for (var i = 0; i < p; i++) {
    var mm = [];
    for (var j = 0; j <= p; j++) {
      mm[ j ] = 0;
    }
    m[ i ] = mm;
  }

  var mpc = [ n ];

  for (var i = 1; i < rs; i++) {
    mpc[ i ] = 0;
  }

  for (i = 0; i < n; i++) {
    var pr = this.matrix[ i ];
    // Process precalculation array
    for (r = 1; r < rs; r++) {
      mpc[ r ] += Math.pow(pr.x, r);
    }
    // Process RH column cells
    m[ 0 ][ p ] += pr.y;
    for (r = 1; r < p; r++) {
      m[ r ][ p ] += Math.pow(pr.x, r) * pr.y;
    }
  }

  // Populate square matrix section
  for (r = 0; r < p; r++) {
    for (c = 0; c < p; c++) {
      m[ r ][ c ] = mpc[ r + c ];
    }
  }

  this.gaussJordanEchelonize(m);

  var terms = [];

  for (var i = 0, len = m.length; i < len; i++) {
    terms[ i ] = m[ i ][ p ];
  }

  return terms;

};

/**
 * Using given degree of fitment, return a function that will calculate the y for a given x
 */

Polyfit.prototype.getPolynomial = function(degree) {

  var degree = parseInt(degree, 10);

  if (isNaN(degree) || degree < 0) throw new Error('Degree must be a positive integer');

  var terms = this.computeCoefficients(degree)
    , eqParts = []
    , len = terms.length;

  eqParts.push(terms[ 0 ]);

  for (var i = 1; i < len; i++) {
    eqParts.push(terms[ i ] + ' * Math.pow(x, ' + i + ')');
  }

  var expr = 'return ' + eqParts.join(' + ') + ';';

  return new Function('x', expr);

};

/**
 * Convert the polynomial to a string expression, mostly useful for visual debugging
 */

Polyfit.prototype.toExpression = function(degree) {

  var degree = parseInt(degree, 10);

  if (isNaN(degree) || degree < 0) throw new Error('Degree must be a positive integer');

  var terms = this.computeCoefficients(degree)
    , eqParts = []
    , len = terms.length;

  eqParts.push(terms[ 0 ]);

  for (var i = 1; i < len; i++) {
    eqParts.push(terms[ i ] + 'x^' + i);
  }

  return eqParts.join(' + ');

};
