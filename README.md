polyfit.js
============

Polynomial fitment solver library for node.js

Heavily, heavily inspired by Paul Lutus (lutusp@arachnoid.com) - http://www.arachnoid.com/polysolve
This was just a way to node-ize / component-ize this library, and format it a little better.

Usage
============

To get a polynomial function that you can pass an x value into
and get the corresponding y value (the parameter is the number of degrees)

```javascript

  var Polyfit = require('polyfit');
  var poly = new Polyfit([ 1, 2, 3, 4, 5 ], [ 0.01, 0.03, -0.02, 0.03, 0.02 ]);
  var solver = poly.getPolynomial(6);

```

Solver will be a function that you can call with an x value.

```javascript

  console.log(solver(1.17));

```

Computing coefficient terms:

```javascript

  var terms = poly.computeCoefficients(6);

```

Use those terms to get a standard error:

```javascript

  var standardError = poly.standardError(terms);

```

Or to get a correlation coefficient

```javascript

  var cc = poly.correlationCoefficient(terms);

```

License
============
(Based on GNU license from Paul Lutus)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.     
                                                                      
You should have received a copy of the GNU General Public License
along with this program; if not, write to the
Free Software Foundation, Inc.,
59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
