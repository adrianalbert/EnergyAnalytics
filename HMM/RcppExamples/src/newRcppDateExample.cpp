// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppParamsExample.h: Rcpp R/C++ interface class library RcppDate example
//
// Copyright (C) 2009 - 2011 Dirk Eddelbuettel and Romain Francois
//
// This file is part of Rcpp.
//
// Rcpp is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Rcpp is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Rcpp.  If not, see <http://www.gnu.org/licenses/>.

#include <Rcpp.h>

RcppExport SEXP newRcppDateExample(SEXP dvsexp, SEXP dtvsexp) {

    try {					// or use BEGIN_RCPP macro

	Rcpp::DateVector dv(dvsexp);
	Rcpp::DatetimeVector dtv(dtvsexp);
	Rcpp::Function formatDate("format.Date");
	Rcpp::Function formatDatetime("format.POSIXct");

	Rprintf("\nIn C++, seeing the following date value\n");
	for (int i=0; i<dv.size(); i++) {
	    Rcpp::Rcout << Rcpp::as<std::string>(formatDate(Rcpp::wrap(dv[i]))) << std::endl;
	    dv[i] = dv[i] + 7;		// shift a week
	}
	Rprintf("\nIn C++, seeing the following datetime value\n");
	for (int i=0; i<dtv.size(); i++) {
	    Rcpp::Rcout << Rcpp::as<std::string>(formatDatetime(Rcpp::wrap(dtv[i]))) << std::endl;
	    dtv[i] = dtv[i] + 0.250;    // shift 250 millisec
	}

	// Build result set to be returned as a list to R.
	return Rcpp::List::create(Rcpp::Named("date",   dv),
				  Rcpp::Named("datetime", dtv));

    } catch( std::exception &ex ) {		// or use END_RCPP macro
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}


