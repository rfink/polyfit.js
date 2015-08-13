"use strict";

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
  if (!(x instanceof Array && y instanceof Array))
  if (!(x instanceof Float64Array && y instanceof Float64Array)) {
      throw new Error('x and y must be arrays');
  }

  // Make sure we have equal lengths
  if (x.length !== y.length) {
    throw new Error('x and y must have the same length');
  }

  this.x = x;
  this.y = y;

}

/**
 * Perform gauss-jordan division
 */

Polyfit.gaussJordanDivide = function gaussJordanDivide(matrix, row, col, numCols) {

  for (var i = col + 1; i < numCols; i++) {
    matrix[ row ][ i ] /= matrix[ row ][ col ];
  }

  matrix[ row ][ col ] = 1;

};

/**
 * Perform gauss-jordan elimination
 */

Polyfit.gaussJordanEliminate = function gaussJordanEliminate(matrix, row, column, numRows, numCols) {

  for (var i = 0; i < numRows; i++) {
    if (i != row && matrix[ i ][ column ] !== 0) {
      for (var j = column + 1; j < numCols; j++) {
        matrix[ i ][ j ] -= matrix[ i ][ column ] * matrix[ row ][ j ];
      }
      matrix[ i ][ column ] = 0;
    }
  }

};

/**
 * Perform gauss-jordan echelon method
 */

Polyfit.gaussJordanEchelonize = function gaussJordanEchelonize(matrix) {

  var rows = matrix.length;
  var cols = matrix[ 0 ].length;
  var i = 0;
  var j = 0;
  var k;
  var swap;

  while (i < rows && j < cols) {
    k = i;
    // Look for non-zero entries in col j at or below row i
    while (k < rows && matrix[ k ][ j ] === 0) {
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
        Polyfit.gaussJordanDivide(matrix, i, j, cols);
      }
      // Eliminate all other non-zero entries
      Polyfit.gaussJordanEliminate(matrix, i, j, rows, cols);
      i++;
    }
    j++;
  }

  return matrix;

};

/**
 * Perform regression
 */

Polyfit.regress = function regress(x, terms) {

  var a = 0;
  var exp = 0;

  for (var i = 0, len = terms.length; i < len; i++) {
    a += terms[ i ] * Math.pow(x, exp++);
  }

  return a;

};

/**
 * Compute correlation coefficient
 */

Polyfit.prototype.correlationCoefficient = function correlationCoefficient(terms) {

  var r = 0;
  var n = this.x.length;
  var sx = 0;
  var sx2 = 0;
  var sy = 0;
  var sy2 = 0;
  var sxy = 0;
  var x;
  var y;

  for (var i = 0; i < n; i++) {
    x = Polyfit.regress(this.x[i], terms);
    y = this.y[i];
    sx += x;
    sy += y;
    sxy += x * y;
    sx2 += x * x;
    sy2 += y * y;
  }

  var div = Math.sqrt((sx2 - (sx * sx) / n) * (sy2 - (sy * sy) / n));

  if (div !== 0) {
    r = Math.pow((sxy - (sx * sy) / n) / div, 2);
  }

  return r;

};

/**
 * Run standard error function
 */

Polyfit.prototype.standardError = function standardError(terms) {

  var r = 0;
  var n = this.x.length;

  if (n > 2) {
    var a = 0;
    for (var i = 0; i < n; i++) {
      a += Math.pow((Polyfit.regress(this.x[i], terms) - this.y[i]), 2);
    }
    r = Math.sqrt(a / (n - 2));
  }

  return r;

};

/**
 * Compute coefficients for given data matrix
 */

Polyfit.prototype.computeCoefficients = function computeCoefficients(p) {

  var n = this.x.length;
  var r;
  var c;
  var rs = 2 * (++p) - 1;
  var i;

  var m = [];

  // Initialize array with 0 values
  for (i = 0; i < p; i++) {
    var mm = [];
    for (var j = 0; j <= p; j++) {
      mm[ j ] = 0;
    }
    m[ i ] = mm;
  }

  var mpc = [ n ];

  for (i = 1; i < rs; i++) {
    mpc[ i ] = 0;
  }

  for (i = 0; i < n; i++) {
    var x = this.x[i];
    var y = this.y[i];

    // Process precalculation array
    for (r = 1; r < rs; r++) {
      mpc[ r ] += Math.pow(x, r);
    }
    // Process RH column cells
    m[ 0 ][ p ] += y;
    for (r = 1; r < p; r++) {
      m[ r ][ p ] += Math.pow(x, r) * y;
    }
  }

  // Populate square matrix section
  for (r = 0; r < p; r++) {
    for (c = 0; c < p; c++) {
      m[ r ][ c ] = mpc[ r + c ];
    }
  }

  Polyfit.gaussJordanEchelonize(m);

  var terms = [];

  for (i = m.length - 1; i >= 0; i--) {
    terms[ i ] = m[ i ][ p ];
  }

  return terms;

};

/**
 * Using given degree of fitment, return a function that will calculate the y for a given x
 */

Polyfit.prototype.getPolynomial = function getPolynomial(degree) {

  degree = parseInt(degree, 10);

  if (isNaN(degree) || degree < 0) throw new Error('Degree must be a positive integer');

  var terms = this.computeCoefficients(degree);
  var eqParts = [];

  eqParts.push(terms[ 0 ]);

  for (var i = 1, len = terms.length; i < len; i++) {
    eqParts.push(terms[ i ] + ' * Math.pow(x, ' + i + ')');
  }

  var expr = 'return ' + eqParts.join(' + ') + ';';

  return new Function('x', expr);

};

/**
 * Convert the polynomial to a string expression, mostly useful for visual debugging
 */

Polyfit.prototype.toExpression = function toExpression(degree) {

  degree = parseInt(degree, 10);

  if (isNaN(degree) || degree < 0) throw new Error('Degree must be a positive integer');

  var terms = this.computeCoefficients(degree);
  var eqParts = [];
  var len = terms.length;

  eqParts.push(terms[ 0 ]);

  for (var i = 1; i < len; i++) {
    eqParts.push(terms[ i ] + 'x^' + i);
  }

  return eqParts.join(' + ');

};
