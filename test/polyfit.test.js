
/* ~ Deps ~ */

var should = require('should');
var Polyfit = require('..');

/* ~ Test data ~ */

var xd = [-1, 0, 1, 2, 3, 5, 7, 9];
var yd = [-1, 3, 2.5, 5, 4, 2, 5, 4];
var degree = 6;

describe('polyfit', function() {

  it('should throw an error if arrays not used', function(done) {
    try {
      var pf = new Polyfit();
    } catch (e) {
      done();
    }
  });
  
  it('should accept Float64Arrays', function(done) {
    var pf = new Polyfit(new Float64Array([1, 0, 1]), new Float64Array([2, 3, 4]));
	done();
  });

  it('should throw an error if arrays not same length', function(done) {
    try {
      var pf = new Polyfit([1], [2, 3]);
    } catch (e) {
      done();
    }
  });

  it('should produce a correct equation', function(done) {
    var pf = new Polyfit(xd, yd);
    var eq = pf.getPolynomial(degree);
    eq(2).should.equal(4.08111112545548);
    eq(3).should.equal(4.502517353251342);
    // TODO: Add assertions
    done();
  });

  it('should produce a valid expression', function(done) {
    var pf = new Polyfit(xd, yd);
    var exp = pf.toExpression(degree);
    exp.should.equal('2.6937037085228717 + 0.9585108884477604x^1 + -1.150528829693737x^2 + 1.0886762123312619x^3 + -0.38856236522551396x^4 + 0.054575046507659646x^5 + -0.002598631007421001x^6');
    done();
  });

  it('should compute correlation coefficient correctly', function(done) {
    var pf = new Polyfit(xd, yd);
    var terms = pf.computeCoefficients(degree);
    pf.correlationCoefficient(terms).should.equal(0.9348988507857894);
    done();
  });

  it('should compute standard error correctly', function(done) {
    var pf = new Polyfit(xd, yd);
    var terms = pf.computeCoefficients(degree);
    pf.standardError(terms).should.equal(0.5434414879841104);
    done();
  });

});