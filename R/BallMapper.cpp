#include <Rcpp.h>
#include <vector>
#include <iostream>
#include <cmath>
using namespace Rcpp;
using namespace std;

/*
//For standard, not symmetric data:
library(Rcpp)
sourceCpp('BallMapper.cpp')

pts <- as.data.frame( read.csv('../circle') )
values = pts[,1]

BallMapperCpp <- function( points , values , epsilon )
{
  output <- SimplifiedBallMapperCppInterface( points , values , epsilon )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp

BallMapperCpp( pts,values,0.3 )
*/

inline double compute_distance
  ( const std::vector< NumericVector >& pts , size_t f , size_t s ,  double p = 2 )
{
  double result = 0;
  for ( size_t i = 0 ; i != pts.size() ; ++i )
  {
     result += pow( ( pts[i][f]-pts[i][s] ) , p );
  }
  return pow(result,1/p);
}//Euclidean_distance
//
//
//
//
//
//
//
//
//
void compute_landmarks
(
 const std::vector< NumericVector >& points ,
 std::vector< std::vector<size_t> >& coverage,
 std::vector< size_t >& landmarks ,
 double epsilon , int number_of_points
)
{
  bool dbg = false;
  int first_not_covered_point = 0;
  int number_of_landmark = -1;
  while ( first_not_covered_point < number_of_points )
  {
     if (dbg) Rcerr << "first_not_covered_point : " << first_not_covered_point << endl;
     landmarks.push_back( first_not_covered_point+1 );//!!!
     ++number_of_landmark;
     //check what is covered by points_vec[first_not_covered_point]:
     for ( int i = 0 ; i != number_of_points ; ++i )
     {
       if ( compute_distance( points , first_not_covered_point , i ) <= epsilon )
       {
         coverage[i].push_back( number_of_landmark+1 );
       }
     }
     //now find the next not covered point:
     while ( (first_not_covered_point!=number_of_points) && (coverage[first_not_covered_point].size() != 0) )
     {
       ++first_not_covered_point;
     }
  }
}
//
//
//
//
//
//
//
//
//
void compute_landmarks_symmetric_point_clouds
(
 const std::vector< NumericVector >& points ,
 std::vector< std::vector<size_t> >& coverage,
 std::vector< size_t >& landmarks ,
 std::vector< std::vector< int > > orbit_of_a_point,
 double epsilon , int number_of_points
)
{
  bool dbg = false;
  int first_not_covered_point = 0;
  int number_of_landmark = -1;
  while ( first_not_covered_point < number_of_points )
  {
     if (dbg) Rcerr << "first_not_covered_point : " << first_not_covered_point << endl;


     //instead of putting a single point to the landmarks list, we will put there the
     //whole orbit of a point. We assume that first_not_covered_point+1 is in its orbit,
     //and this is why we do not need the two lines below:
     //landmarks.push_back( first_not_covered_point+1 );//!!!
     //++number_of_landmark;
     for ( size_t j = 0 ; j != orbit_of_a_point[ first_not_covered_point+1 ].size() ; ++j )
     {
         landmarks.push_back( orbit_of_a_point[ first_not_covered_point ][j] );
         ++number_of_landmark;

         //check what is covered by points_vec[first_not_covered_point]:
         for ( int i = 0 ; i != number_of_points ; ++i )
         {
             if ( compute_distance( points , orbit_of_a_point[ first_not_covered_point ][j] , i ) <= epsilon )
             {
                coverage[i].push_back( number_of_landmark+1 );
             }
         }
     }

     //now find the next not covered point:
     while ( (first_not_covered_point!=number_of_points) && (coverage[first_not_covered_point].size() != 0) )
     {
       ++first_not_covered_point;
     }
  }

  //now we still should remove repetitions from coverage list -- it may show up here, as we add multiple landmarks
  //at the same time.
  for ( size_t i = 0 ; i != coverage.size() ; ++i )
  {
      //remove repetitions from coverage[i]:
      std::sort( coverage[i].begin() , coverage[i].end() );
      coverage[i].erase( std::unique( coverage[i].begin() , coverage[i].end() ) , coverage[i].end() );
  }
}













//internal procedure
void internal_procedure_fill_points_covered_by_landmarks
( 
NumericVector& numer_of_covered_points , 
std::vector< std::vector< int > >& points_covered_by_landmarks , 
const std::vector< size_t >& landmarks , 
const std::vector< std::vector<size_t> >& coverage 
)
{
	bool dbg = false;
	//Now we will compute points_covered_by_landmarks. Firstly let us initialize all the structures:
	numer_of_covered_points = NumericVector( landmarks.size() , 2 );//We initialize it to 2, as othervise we get very unbalanced balls sizes.
	for ( size_t i = 0 ; i != coverage.size() ; ++i )
	{
		for ( size_t j = 0 ; j != coverage[i].size() ; ++j )
		{
			numer_of_covered_points[ coverage[i][j]-1 ]++;
		}
	}

	points_covered_by_landmarks = std::vector< std::vector< int > > ( landmarks.size() );
	for ( size_t i = 0 ; i != landmarks.size() ; ++i )
	{
		std::vector<int> aa;
		aa.reserve( numer_of_covered_points[i] );
		points_covered_by_landmarks[i] = aa;
	}
	//now when the variables are initialized, we can fill in the numer_of_covered_points list:
	for ( size_t i = 0 ; i != coverage.size() ; ++i )
	{
		for ( size_t j = 0 ; j != coverage[i].size() ; ++j )
		{
			points_covered_by_landmarks[ coverage[i][j]-1 ].push_back( i+1 );
		}
	}
	if (dbg) Rcerr << "points_covered_by_landmarks.size() : " << points_covered_by_landmarks.size() << endl;
}//internal_procedure_fill_points_covered_by_landmarks
//
//
//
//
//
//
//
//
//
void
internal_procedure_fill_coloring
(
std::vector< double >& coloring , 
const std::vector< std::vector< int > >& points_covered_by_landmarks , 
NumericVector& values
)
{
	double dbg = false;
	coloring = std::vector< double >( points_covered_by_landmarks.size() , 0 );
	for ( size_t i = 0 ; i != points_covered_by_landmarks.size() ; ++i )
	{
		double av = 0;
		for ( size_t j = 0 ; j != points_covered_by_landmarks[i].size() ; ++j )
		{
			av += values[ points_covered_by_landmarks[i][j]-1 ];
		}
		av = av / (double)points_covered_by_landmarks[i].size();
		coloring[i] = av;
	}

	if (dbg)
	{
		Rcerr << "Here is the coloring : \n";
		for ( size_t i = 0 ; i != coloring.size() ; ++i )
		{
			Rcerr << coloring[i] << " , ";
		}
	}
}//internal_procedure_fill_coloring
//
//
//
//
//
//
//
//
//
void
internal_procedure_build_graph
(
std::vector< std::vector< int > >& graph_incidence_matrix,
NumericVector& from,
NumericVector& to,
NumericVector& strength_of_edges,
const std::vector< size_t >& landmarks,
const std::vector< std::vector<size_t> > coverage
)
{
  bool dbg = false;
  graph_incidence_matrix = std::vector< std::vector< int > >( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
    graph_incidence_matrix[i] = std::vector<int>( i );
  }


  if (dbg) Rcerr << "coverage.size() : " << coverage.size() << endl;

  for ( size_t i = 0 ; i < coverage.size() ; ++i )
  {
    for ( size_t j = 0 ; j < coverage[i].size() ; ++j )
    {
      for ( size_t k = j+1 ; k < coverage[i].size() ; ++k )
      {
        //note that the landmarks in coverage are sorted from smallest to largest
        //therefore we only need this:
        //Rcerr << coverage[i][k] << " " << coverage[i][j] << endl;
        graph_incidence_matrix[ coverage[i][k]-1 ][ coverage[i][j]-1  ]++;
      }
    }
  }


  //first let us count the number of edges in the graph:
  int number_of_edges = 0;
  for ( size_t i = 0 ; i != graph_incidence_matrix.size() ; ++i )
  {
    for ( size_t j = 0 ; j != graph_incidence_matrix[i].size() ; ++j )
    {
       if ( graph_incidence_matrix[i][j] != 0 )++number_of_edges;
    }
  }
  if (dbg)Rcerr << "Number of edges in the graph : " << number_of_edges << endl;

  from = NumericVector(number_of_edges);
  to = NumericVector(number_of_edges);
  strength_of_edges = NumericVector(number_of_edges);

  int edg_no = 0;

  for ( size_t i = 0 ; i != graph_incidence_matrix.size() ; ++i )
  {
    for ( size_t j = 0 ; j != graph_incidence_matrix[i].size() ; ++j )
    {
      if ( graph_incidence_matrix[i][j] != 0 )
      {
         from[edg_no] = j+1;
         to[edg_no] = i+1;
         strength_of_edges[edg_no] = graph_incidence_matrix[i][j];
         ++edg_no;
      }
    }
  }
}//internal_procedure_build_graph
//
//
//
//
//
//
//
//
//
/*
//For standard, not symmetric data:
library(Rcpp)
library(BallMapper)

sourceCpp('BallMapper.cpp')

pts <- as.data.frame( read.csv('circle', header=FALSE)  )
values = pts[,1]

BallMapperCpp <- function( points , values , epsilon )
{
  output <- BallMapperCppInterface( points , values , epsilon )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp

bm <- BallMapperCpp( pts,values,0.4 )
BallMapper::ColorIgraphPlot(bm)
*/
// [[Rcpp::export]]
List BallMapperCppInterface( const DataFrame& points_df, const DataFrame& values_df , double epsilon  )
{
  bool dbg = false;
  if ( points_df.size() == 0 )
  {
    cerr << "No points in the BallMapperCpp procedure, the program will now terminate";
    throw "No points in the BallMapperCpp procedure, the program will now terminate";
  }

  std::vector< NumericVector > points( points_df.size() );
  for ( int i = 0 ; i != points_df.size() ; ++i )
  {
    points[i] = points_df[i];
  }
  int number_of_points = points[0].size();
  if (dbg) Rcerr << "Number of points : " << number_of_points << endl;

  NumericVector values = values_df[0];
  if (dbg) Rcerr << "values.size() : " << values.size() << endl;

  std::vector< std::vector<size_t> > coverage( number_of_points );
  std::vector< size_t > landmarks;
  landmarks.reserve( (size_t)(0.2*number_of_points) );

  //here we outsource computations of landmark points:
  compute_landmarks( points , coverage, landmarks , epsilon , number_of_points );

  if (dbg) Rcerr << "landmarks.size() : " << landmarks.size() << endl;
  if (dbg) Rcerr << "coverage.size() : " << coverage.size() << endl;

  if (dbg)
  {
    Rcerr << "Here are the landmarks: \n";
    for ( size_t i = 0 ; i != landmarks.size() ; ++i )
    {
      Rcerr << landmarks[i] << " , ";
    }
  }

//change1
  NumericVector numer_of_covered_points;
  std::vector< std::vector< int > > points_covered_by_landmarks;
  internal_procedure_fill_points_covered_by_landmarks( numer_of_covered_points , points_covered_by_landmarks ,  landmarks , coverage );

  //now let us deal with the coloring of the vertices:
  //change2
  std::vector< double > coloring;
  internal_procedure_fill_coloring( coloring , points_covered_by_landmarks , values ); 
 
 
  //Now let us create the graph and record the strength of each edge. To measure this, we will create the incidence matrix of the graph.
  //change3
  std::vector< std::vector< int > > graph_incidence_matrix;
  NumericVector from;
  NumericVector to;
  NumericVector strength_of_edges;
  internal_procedure_build_graph( graph_incidence_matrix, from, to, strength_of_edges, landmarks, coverage );
  
  NumericVector verts( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
     verts[i] = i+1;
  }


  List ret;
  ret["vertices"] = cbind(verts,numer_of_covered_points);
  ret["edges"] = cbind(from,to);
  ret["strength_of_edges"] = strength_of_edges;
  ret["points_covered_by_landmarks"] = points_covered_by_landmarks;
  ret["landmarks"] = landmarks;
  ret["coloring"] = coloring;
  ret["coverage"] = coverage;
  //ret["numer_of_covered_points"] = numer_of_covered_points;
  return ret;
}

//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
inline double compute_distance_standard_points
  ( const std::vector< double >& pt1 , const std::vector< double >& pt2 , double p = 2 )
{
  double result = 0;
  for ( size_t i = 0 ; i != pt1.size() ; ++i )
  {
     result += pow( ( pt1[i]-pt2[i] ) , p );
  }
  return pow(result,1/p);
}//compute_distance_standard_points
//
//
//
//
//
//
//
//
void compute_landmarks_not_transposed_pts_group_action
                                         ( std::vector< std::vector<size_t> >& coverage ,
                                           std::vector< size_t > & landmarks ,
                                           const std::vector< std::vector<double> >& points ,
                                           const std::vector< std::vector<int> >& orbit ,
                                           double epsilon
                                         )
{
   bool dbg = false;

   //here we outsource computations of landmark points:
  size_t current_point = 0;
  size_t current_landmark = 0;

  if ( dbg )
  {
     Rcerr << "orbit.size() : " << orbit.size() << std::endl;
  }

  while ( true )
  {
      while ( (current_point != points.size()) && (coverage[current_point].size() != 0) )
      {
          ++current_point;
      }

      if ( dbg )Rcerr << "Current point : " << current_point << std::endl;

      if ( current_point == points.size() )break;
      for ( size_t i = 0 ; i != orbit[ current_point ].size() ; ++i )
      {
           if ( dbg )Rcerr << "orbit["<<current_point<<"][" << i << "] : " << orbit[ current_point ][i] << endl;

           landmarks.push_back( orbit[ current_point ][i]-1 );
           for ( size_t j = 0 ; j != points.size() ; ++j )
           {
              if ( compute_distance_standard_points( points[j] , points[ orbit[ current_point ][i]-1 ] ) <= epsilon )
              {
                  coverage[j].push_back( current_landmark+1 );
              }
           }
           current_landmark++;
      }
      if ( dbg )Rcerr << "Out of the internal while loop. \n";
  }
}
//
//
//
//
//
//
//
//
//
void compute_landmarks_not_transposed_pts( std::vector< std::vector<size_t> >& coverage ,
                                           std::vector< size_t > & landmarks ,
                                           const std::vector< std::vector<double> >& points ,
                                           double epsilon
                                         )
{
   //here we outsource computations of landmark points:
  size_t current_point = 0;
  size_t current_landmark = 0;
  while ( true )
  {
      while ( (current_point != points.size()) && (coverage[current_point].size() != 0) )
      {
          ++current_point;
      }
      if ( current_point == points.size() )break;
      landmarks.push_back( current_point );
      for ( size_t j = 0 ; j != points.size() ; ++j )
      {
          if ( compute_distance_standard_points( points[j] , points[current_point] ) <= epsilon )
          {
              coverage[j].push_back( current_landmark+1 );
          }
      }
      current_landmark++;
  }
}

std::vector< std::vector<double> > transpose_points_from_R( const DataFrame& points_df )
{
    std::vector< NumericVector > points_transposed( points_df.size() );
    for ( int i = 0 ; i != points_df.size() ; ++i )
    {
        points_transposed[i] = points_df[i];
    }

    std::vector< std::vector<double> > points( points_transposed[0].size() );

    size_t dim_of_points = points_transposed.size();
    for ( int i = 0 ; i != points_transposed[0].size() ; ++i )
    {
        std::vector< double > pt( dim_of_points );
        for ( size_t j = 0 ; j != dim_of_points ; ++j )
        {
         pt[ j ] = points_transposed[j][i];
        }
        points[i] = pt;
    }
    return points;
}


std::vector< std::vector<int> > transpose_points_from_R_int_version( const DataFrame& points_df )
{
    std::vector< NumericVector > points_transposed( points_df.size() );
    for ( int i = 0 ; i != points_df.size() ; ++i )
    {
        points_transposed[i] = points_df[i];
    }

    std::vector< std::vector<int> > points( points_transposed[0].size() );

    size_t dim_of_points = points_transposed.size();
    for ( int i = 0 ; i != points_transposed[0].size() ; ++i )
    {
        std::vector< int > pt( dim_of_points );
        for ( size_t j = 0 ; j != dim_of_points ; ++j )
        {
         pt[ j ] = points_transposed[j][i];
        }
        points[i] = pt;
    }
    return points;
}


/*
library(Rcpp)
library(BallMapper)

sourceCpp('BallMapper.cpp')

pts <- as.data.frame( read.csv('circle', header=FALSE)  )
values = pts[,1]

BallMapperCpp <- function( points , values , epsilon )
{
  output <- SimplifiedBallMapperCppInterface( points , values , epsilon )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp

bm <- BallMapperCpp( pts,values,0.4 )
BallMapper::ColorIgraphPlot(bm)
*/
// [[Rcpp::export]]
List SimplifiedBallMapperCppInterface( const DataFrame& points_df , const DataFrame& values_df , double epsilon  )
{
  bool dbg = false;
  if ( points_df.size() == 0 )
  {
    Rcerr << "No points in the BallMapperCpp procedure, the program will now terminate";
    throw "No points in the BallMapperCpp procedure, the program will now terminate";
  }

  //the points we obtain from R are unnaturally transposed, so we transpose them back to the situation when we
  //have one point per row.
  //First we need to store them as vector of NumericVectorS:
  std::vector< std::vector<double> > points = transpose_points_from_R( points_df );

  int number_of_points = points.size();
  if (dbg) Rcerr << "Number of points : " << number_of_points << endl;

  NumericVector values = values_df[0];
  if (dbg) Rcerr << "values.size() : " << values.size() << endl;


  std::vector< std::vector<size_t> > coverage( number_of_points );
  std::vector< size_t > landmarks;
  landmarks.reserve( (size_t)(0.2*number_of_points) );

  compute_landmarks_not_transposed_pts( coverage , landmarks , points , epsilon );


  if (dbg) Rcerr << "landmarks.size() : " << landmarks.size() << endl;
  if (dbg) Rcerr << "coverage.size() : " << coverage.size() << endl;

  if (dbg)
  {
    Rcerr << "Here are the landmarks: \n";
    for ( size_t i = 0 ; i != landmarks.size() ; ++i )
    {
      Rcerr << landmarks[i] << " , ";
    }
  }


  NumericVector numer_of_covered_points;
  std::vector< std::vector< int > > points_covered_by_landmarks;
  internal_procedure_fill_points_covered_by_landmarks( numer_of_covered_points , points_covered_by_landmarks , landmarks , coverage );

  //now let us deal with the coloring of the vertices:
  std::vector< double > coloring;
  internal_procedure_fill_coloring( coloring , points_covered_by_landmarks , values ); 

  std::vector< std::vector< int > > graph_incidence_matrix;
  NumericVector from;
  NumericVector to;
  NumericVector strength_of_edges;
  internal_procedure_build_graph( graph_incidence_matrix, from, to, strength_of_edges, landmarks, coverage );

  NumericVector verts( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
     verts[i] = i+1;
  }

  List ret;
  ret["vertices"] = cbind(verts,numer_of_covered_points);
  ret["edges"] = cbind(from,to);
  ret["strength_of_edges"] = strength_of_edges;
  ret["points_covered_by_landmarks"] = points_covered_by_landmarks;
  ret["landmarks"] = landmarks;
  ret["coloring"] = coloring;
  ret["coverage"] = coverage;
  return ret;

return 1;
}



/*
To get the landmark points in R:
aa <- compute_landmark_indices(pts,0.4)
aa = t(aa)
lands <- pts[aa,]
and then they can be used to call the next subroutine:
*/
//[[Rcpp::export]]
NumericVector compute_landmark_indices
(
 const std::vector< NumericVector >& points ,
 double epsilon
)
{
  int number_of_points = points[0].size();
  bool dbg = false;
  int first_not_covered_point = 0;
  int number_of_landmark = -1;

  //let us initialize the landmarks:
  NumericVector landmarks;

  //and coverage:
  std::vector< std::vector<size_t> > coverage( number_of_points );

  while ( first_not_covered_point < number_of_points )
  {
     if (dbg) Rcerr << "first_not_covered_point : " << first_not_covered_point << endl;
     landmarks.push_back( first_not_covered_point+1 );//!!!
     ++number_of_landmark;
     //check what is covered by points_vec[first_not_covered_point]:
     for ( int i = 0 ; i != number_of_points ; ++i )
     {
       if ( compute_distance( points , first_not_covered_point , i ) <= epsilon )
       {
         coverage[i].push_back( number_of_landmark+1 );
       }
     }
     //now find the next not covered point:
     while ( (first_not_covered_point!=number_of_points) && (coverage[first_not_covered_point].size() != 0) )
     {
       ++first_not_covered_point;
     }
  }
  return landmarks;
}

/*
new_epsilon <- computeDistanceBetweenPtClouds(pts,lands)
*/
// [[Rcpp::export]]
double computeDistanceBetweenPtClouds
(
  const DataFrame& pts1_df ,
  const DataFrame& pts2_df
)
{
  std::vector< std::vector<double> > pts1 = transpose_points_from_R( pts1_df );
  std::vector< std::vector<double> > pts2 = transpose_points_from_R( pts2_df );
  double max_distance = 0;
  double d;
  //for every point in one point cloud:
  for ( size_t i = 0 ; i != pts1.size() ; ++i )
  {
      //compute the minimal distance to points from pts2:
      double min_distance = std::numeric_limits<double>::infinity();
      for ( size_t j = 0 ; j != pts2.size() ; ++j )
      {
          d = compute_distance_standard_points( pts1[i] , pts2[j] );
          if ( d < min_distance )min_distance = d;
      }
      if ( min_distance > max_distance )max_distance = min_distance;
  }

  //and now the symmetric version:

  //for every point in second point cloud:
  //for ( size_t i = 0 ; i != pts2.size() ; ++i )
  //{
  //    //compute the minimal distance to points from pts1:
  //    double min_distance = std::numeric_limits<double>::infinity();
  //    for ( size_t j = 0 ; j != pts1.size() ; ++j )
  //    {
  //        d = compute_distance_standard_points( pts2[i] , pts1[j] );
  //        if ( d < min_distance )min_distance = d;
  //    }
  //    if ( min_distance > max_distance )max_distance = min_distance;
  //}
  return max_distance;
}//computeDistanceBetweenPtClouds


/*
library(Rcpp)
sourceCpp('BallMapper.cpp')

pts <- as.data.frame( read.csv('../circle') )
values = pts[,1]

aa <- compute_landmark_indices(pts,0.4)
aa = t(aa)
lands <- pts[aa,]

new_epsilon <- computeDistanceBetweenPtClouds(pts,lands)

BallMapperLandCpp <- function( points , landmarks , values , epsilon )
{
  output <- SimplifiedBallMapperCppInterfaceBasedOnLands( points , landmarks , values , epsilon )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperLandCpp

BallMapperLandCpp( pts,lands,values,0.4 )
*/
/*
library(Rcpp)
library(BallMapper)

sourceCpp('BallMapper.cpp')

pts <- as.data.frame( read.csv('circle', header=FALSE)  )
values = pts[,1]

BallMapperCpp <- function( points , values , epsilon )
{
  aa <- compute_landmark_indices(pts,0.4)
  aa = t(aa)
  lands <- pts[aa,]
  output <- SimplifiedBallMapperCppInterfaceBasedOnLands( points , values , lands , epsilon )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp

bm <- BallMapperCpp( pts,values,0.4 )
BallMapper::ColorIgraphPlot(bm)
*/
// [[Rcpp::export]]
List SimplifiedBallMapperCppInterfaceBasedOnLands
(
  const DataFrame& points_df ,
  const DataFrame& landmarks_df ,
  const DataFrame& values_df ,
  double epsilon
)
{
  bool dbg = false;
  if ( points_df.size() == 0 )
  {
    Rcerr << "No points in the BallMapperCpp procedure, the program will now terminate";
    throw "No points in the BallMapperCpp procedure, the program will now terminate";
  }

  std::vector< std::vector<double> > points = transpose_points_from_R( points_df );

  int number_of_points = points.size();
  if (dbg) Rcerr << "Number of points : " << number_of_points << endl;

  NumericVector values = values_df[0];
  if (dbg) Rcerr << "values.size() : " << values.size() << endl;

  std::vector< std::vector<double> > landmarks = transpose_points_from_R( landmarks_df );
  if (dbg) Rcerr << "landmarks.size() : " << landmarks.size() << endl;

  std::vector< std::vector<size_t> > coverage( number_of_points );
  for ( size_t i = 0 ; i != points.size() ; ++i )
  {
      for ( size_t j = 0 ; j != landmarks.size() ; ++j )
      {
          if ( compute_distance_standard_points( points[i] , landmarks[j] ) <= epsilon )
          {
              coverage[i].push_back( j+1 );
          }
      }
  }
  

  //Now we will compute points_covered_by_landmarks. Firstly let us initialize all the structures:
  NumericVector numer_of_covered_points( landmarks.size() , 2 );//We initialize it to 2, as othervise we get very unbalanced balls sizes.
  for ( size_t i = 0 ; i != coverage.size() ; ++i )
  {
    for ( size_t j = 0 ; j != coverage[i].size() ; ++j )
    {
        numer_of_covered_points[ coverage[i][j]-1 ]++;
    }
    //Rcerr << endl;
  }

  std::vector< std::vector< int > > points_covered_by_landmarks( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
    std::vector<int> aa;
    aa.reserve( numer_of_covered_points[i] );
    points_covered_by_landmarks[i] = aa;
  }

  //now when the variables are initialized, we can fill in the numer_of_covered_points list:
  for ( size_t i = 0 ; i != coverage.size() ; ++i )
  {
    for ( size_t j = 0 ; j != coverage[i].size() ; ++j )
    {
        points_covered_by_landmarks[ coverage[i][j]-1 ].push_back( i+1 );
    }
  }
  if (dbg) Rcerr << "points_covered_by_landmarks.size() : " << points_covered_by_landmarks.size() << endl;


  std::vector< double > coloring;
  internal_procedure_fill_coloring( coloring , points_covered_by_landmarks , values ); 

  //Now let us create the graph and record the strength of each edge. To measure this, we will create the incidence matrix of the graph.
  std::vector< std::vector< int > > graph_incidence_matrix( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
    graph_incidence_matrix[i] = std::vector<int>( i );
  }


  if (dbg) Rcerr << "coverage.size() : " << coverage.size() << endl;

  for ( size_t i = 0 ; i < coverage.size() ; ++i )
  {
    for ( size_t j = 0 ; j < coverage[i].size() ; ++j )
    {
      for ( size_t k = j+1 ; k < coverage[i].size() ; ++k )
      {
        //note that the landmarks in coverage are sorted from smallest to largest
        //therefore we only need this:
        //Rcerr << coverage[i][k] << " " << coverage[i][j] << endl;
        graph_incidence_matrix[ coverage[i][k]-1 ][ coverage[i][j]-1  ]++;
      }
    }
  }


  //first let us count the number of edges in the graph:
  int number_of_edges = 0;
  for ( size_t i = 0 ; i != graph_incidence_matrix.size() ; ++i )
  {
    for ( size_t j = 0 ; j != graph_incidence_matrix[i].size() ; ++j )
    {
       if ( graph_incidence_matrix[i][j] != 0 )++number_of_edges;
    }
  }
  if (dbg)Rcerr << "Number of edges in the graph : " << number_of_edges << endl;

  NumericVector from(number_of_edges);
  NumericVector to(number_of_edges);
  NumericVector strength_of_edges(number_of_edges);

  int edg_no = 0;

  for ( size_t i = 0 ; i != graph_incidence_matrix.size() ; ++i )
  {
    for ( size_t j = 0 ; j != graph_incidence_matrix[i].size() ; ++j )
    {
      if ( graph_incidence_matrix[i][j] != 0 )
      {
         from[edg_no] = j+1;
         to[edg_no] = i+1;
         strength_of_edges[edg_no] = graph_incidence_matrix[i][j];
         ++edg_no;

      }
    }
  }
  

  NumericVector verts( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
     verts[i] = i+1;
  }

  List ret;
  ret["vertices"] = cbind(verts,numer_of_covered_points);
  ret["edges"] = cbind(from,to);
  ret["strength_of_edges"] = strength_of_edges;
  ret["points_covered_by_landmarks"] = points_covered_by_landmarks;
  ret["landmarks"] = landmarks;
  ret["coloring"] = coloring;
  ret["coverage"] = coverage;
  return ret;

return 1;
}//SimplifiedBallMapperCppInterfaceBasedOnLands








/*
//For data with some group action:

library(Rcpp)
sourceCpp('BallMapper.cpp')

points <- read.csv('three_shifted_bloobs.csv',header=FALSE)
values <- points[,1]
orbit <- read.csv('orbit_three_shifted_bloobs.csv',header=FALSE)
epsilon <- 0.2

BallMapperGroupActionCpp <- function( points , values , epsilon , orbit )
{
  output <- SimplifiedBallMapperCppInterfaceGroupAction( points , values , epsilon , orbit )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp

bm <- BallMapperGroupActionCpp( points , values , epsilon , orbit )
BallMapper::ColorIgraphPlot(bm)
*/

// [[Rcpp::export]]
List SimplifiedBallMapperCppInterfaceGroupAction( const DataFrame& points_df , const DataFrame& values_df , double epsilon , const DataFrame& orbit_ )
{
  bool dbg = false;
  if ( points_df.size() == 0 )
  {
    Rcerr << "No points in the BallMapperCpp procedure, the program will now terminate";
    throw "No points in the BallMapperCpp procedure, the program will now terminate";
  }

  //the points we obtain from R are unnaturally transposed, so we transpose them back to the situation when we
  //have one point per row.
  //First we need to store them as vector of NumericVectorS:

  std::vector< std::vector<double> > points = transpose_points_from_R( points_df );

  int number_of_points = points.size();
  if (dbg) Rcerr << "Number of points : " << number_of_points << endl;

  std::vector< std::vector<int> > orbit = transpose_points_from_R_int_version( orbit_ );
  if (dbg) Rcerr << "orbit.size() : " << orbit.size() << endl;


  NumericVector values = values_df[0];
  if (dbg) Rcerr << "values.size() : " << values.size() << endl;


  std::vector< std::vector<size_t> > coverage( number_of_points );
  std::vector< size_t > landmarks;
  landmarks.reserve( (size_t)(0.2*number_of_points) );


  if (dbg) Rcerr << "Entering compute_landmarks_not_transposed_pts_group_action.\n";
  compute_landmarks_not_transposed_pts_group_action( coverage , landmarks ,  points ,  orbit , epsilon );


  if (dbg) Rcerr << "landmarks.size() : " << landmarks.size() << endl;
  if (dbg) Rcerr << "coverage.size() : " << coverage.size() << endl;

  if (dbg)
  {
    Rcerr << "Here are the landmarks: \n";
    for ( size_t i = 0 ; i != landmarks.size() ; ++i )
    {
      Rcerr << landmarks[i] << " , ";
    }
  }

  NumericVector numer_of_covered_points;
  std::vector< std::vector< int > > points_covered_by_landmarks;
  internal_procedure_fill_points_covered_by_landmarks( numer_of_covered_points , points_covered_by_landmarks ,  landmarks , coverage );

  std::vector< double > coloring;
  internal_procedure_fill_coloring( coloring , points_covered_by_landmarks , values ); 

  std::vector< std::vector< int > > graph_incidence_matrix;
  NumericVector from;
  NumericVector to;
  NumericVector strength_of_edges;
  internal_procedure_build_graph( graph_incidence_matrix, from, to, strength_of_edges, landmarks, coverage );


  NumericVector verts( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
     verts[i] = i+1;
  }

  List ret;
  ret["vertices"] = cbind(verts,numer_of_covered_points);
  ret["edges"] = cbind(from,to);
  ret["strength_of_edges"] = strength_of_edges;
  ret["points_covered_by_landmarks"] = points_covered_by_landmarks;
  ret["landmarks"] = landmarks;
  ret["coloring"] = coloring;
  ret["coverage"] = coverage;
  return ret;

return 1;
}


















































































//SPARSE REPRESENTATION CODE BELOW:
std::vector< std::vector< std::pair<unsigned , double > > > transpose_points_from_R_to_sparse_vector( const DataFrame& points_df )
{
    std::vector< NumericVector > points_transposed( points_df.size() );
    for ( int i = 0 ; i != points_df.size() ; ++i )
    {
        points_transposed[i] = points_df[i];
    }

    std::vector< std::vector< std::pair<unsigned , double > > > points( points_transposed[0].size() );

    size_t dim_of_points = points_transposed.size();
    for ( int i = 0 ; i != points_transposed[0].size() ; ++i )
    {
        std::vector< std::pair<unsigned , double > > pt;
        pt.reserve( dim_of_points );
        for ( size_t j = 0 ; j != dim_of_points ; ++j )
        {
        	if ( points_transposed[j][i] != 0 )
        	{        
      		   pt.push_back( std::make_pair( j , points_transposed[j][i] ) );
      		}
        }
        points[i] = pt;
    }
    return points;
}



inline double compute_distance_standard_points_sparse_points
  ( const std::vector< std::pair<unsigned,double> >& pt1 , const std::vector< std::pair<unsigned,double> >& pt2 , double p = 2 )
{
  double result = 0;
  size_t pt1_it = 0;
  size_t pt2_it = 0;
  
  while ( ( pt1_it < pt1.size() ) && ( pt2_it < pt2.size() ) )
  {
  	if ( pt1[pt1_it].first < pt2[pt2_it].first )
  	{
	        //we move pt1:
	        result += pow( ( pt1[pt1_it].second ) , p );
		++pt1_it;  	  		  		
	}
	else
	{
		if ( pt1[pt1_it].first > pt2[pt2_it].first )
		{
			//we move pt1:
	        	result += pow( ( pt2[pt2_it].second ) , p );
			++pt2_it;  	
		}
		else
		{
			//in this case pt1[pt1_it].first == pt2[pt2_it]
			result += pow( ( pt1[pt1_it].second - pt2[pt2_it].second ) , p );			
			++pt1_it;
			++pt2_it;
		}
	}
  }
  
  //if there is anything left, we count it here:
  while ( pt1_it < pt1.size() )
  {
	  result += pow( ( pt1[pt1_it].second ) , p );
	  ++pt1_it;
  }
  while ( pt2_it < pt2.size() )
  {
	result += pow( ( pt2[pt2_it].second ) , p );
	pt2_it++;
  }
  
  return pow(result,1/p);
}//compute_distance_standard_points_sparse_points



void compute_landmarks_not_transposed_pts_group_action_sparse_points
                                         ( std::vector< std::vector<size_t> >& coverage ,
                                           std::vector< size_t > & landmarks ,
                                           const std::vector< std::vector< std::pair<unsigned,double> > >& points ,
                                           const std::vector< std::vector<int> >& orbit ,
                                           double epsilon
                                         )
{
   bool dbg = false;

   //here we outsource computations of landmark points:
  size_t current_point = 0;
  size_t current_landmark = 0;

  if ( dbg )
  {
     Rcerr << "orbit.size() : " << orbit.size() << std::endl;
  }

  while ( true )
  {
      while ( (current_point != points.size()) && (coverage[current_point].size() != 0) )
      {
          ++current_point;
      }

      if ( dbg )Rcerr << "Current point : " << current_point << std::endl;

      if ( current_point == points.size() )break;
      for ( size_t i = 0 ; i != orbit[ current_point ].size() ; ++i )
      {
           if ( dbg )Rcerr << "orbit["<<current_point<<"][" << i << "] : " << orbit[ current_point ][i] << endl;

  
           landmarks.push_back( orbit[ current_point ][i]-1 );
           for ( size_t j = 0 ; j != points.size() ; ++j )
           {
              if ( compute_distance_standard_points_sparse_points( points[j] , points[ orbit[ current_point ][i]-1 ] ) <= epsilon )
              {
                  coverage[j].push_back( current_landmark+1 );
              }
           }
           current_landmark++;
      }
      if ( dbg )Rcerr << "Out of the internal while loop. \n";
  }
}





/*
//For data with some group action AND SPARSE REPRESENTATION:

library(Rcpp)
sourceCpp('BallMapper.cpp')

points <- read.csv('three_shifted_bloobs_HD.csv',header=FALSE)
values <- points[,5]
orbit <- read.csv('orbit_three_shifted_bloobs.csv',header=FALSE)
epsilon <- 0.2

BallMapperGroupActionCpp <- function( points , values , epsilon , orbit )
{
  output <- SimplifiedBallMapperCppInterfaceGroupActionAndSparseRepresentation( points , values , epsilon , orbit )
  colnames(output$vertices) = c('id','size')
  return_list <- output
}#BallMapperCpp

bm <- BallMapperGroupActionCpp( points , values , epsilon , orbit )
BallMapper::ColorIgraphPlot(bm)
*/

// [[Rcpp::export]]
List SimplifiedBallMapperCppInterfaceGroupActionAndSparseRepresentation( const DataFrame& points_df , const DataFrame& values_df , double epsilon , const DataFrame& orbit_ )
{
  bool dbg = false                                ;
  if ( points_df.size() == 0 )
  {
    Rcerr << "No points in the BallMapperCpp procedure, the program will now terminate";
    throw "No points in the BallMapperCpp procedure, the program will now terminate";
  }

  //the points we obtain from R are unnaturally transposed, so we transpose them back to the situation when we
  //have one point per row.
  //First we need to store them as vector of NumericVectorS:

  std::vector< std::vector< std::pair<unsigned,double> > > points = transpose_points_from_R_to_sparse_vector( points_df );
  
//  Rcerr << "points[0].size()  : " << points[0].size() << endl;

  int number_of_points = points.size();
  if (dbg) Rcerr << "Number of points : " << number_of_points << endl;

  std::vector< std::vector<int> > orbit = transpose_points_from_R_int_version( orbit_ );
  if (dbg) Rcerr << "orbit.size() : " << orbit.size() << endl;


  NumericVector values = values_df[0];
  if (dbg) Rcerr << "values.size() : " << values.size() << endl;


  std::vector< std::vector<size_t> > coverage( number_of_points );
  std::vector< size_t > landmarks;
  landmarks.reserve( (size_t)(0.2*number_of_points) );


  if (dbg) Rcerr << "Entering compute_landmarks_not_transposed_pts_group_action.\n";
  compute_landmarks_not_transposed_pts_group_action_sparse_points( coverage , landmarks ,  points ,  orbit , epsilon );


  if (dbg) Rcerr << "landmarks.size() : " << landmarks.size() << endl;
  if (dbg) Rcerr << "coverage.size() : " << coverage.size() << endl;

  if (dbg)
  {
    Rcerr << "Here are the landmarks: \n";
    for ( size_t i = 0 ; i != landmarks.size() ; ++i )
    {
      Rcerr << landmarks[i] << " , ";
    }
  }

  NumericVector numer_of_covered_points;
  std::vector< std::vector< int > > points_covered_by_landmarks;
  internal_procedure_fill_points_covered_by_landmarks( numer_of_covered_points , points_covered_by_landmarks ,  landmarks , coverage );

  std::vector< double > coloring;
  internal_procedure_fill_coloring( coloring , points_covered_by_landmarks , values ); 

  std::vector< std::vector< int > > graph_incidence_matrix;
  NumericVector from;
  NumericVector to;
  NumericVector strength_of_edges;
  internal_procedure_build_graph( graph_incidence_matrix, from, to, strength_of_edges, landmarks, coverage );


  NumericVector verts( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
     verts[i] = i+1;
  }

  List ret;
  ret["vertices"] = cbind(verts,numer_of_covered_points);
  ret["edges"] = cbind(from,to);
  ret["strength_of_edges"] = strength_of_edges;
  ret["points_covered_by_landmarks"] = points_covered_by_landmarks;
  ret["landmarks"] = landmarks;
  ret["coloring"] = coloring;
  ret["coverage"] = coverage;
  return ret;

return 1;
}






/*
 library(Rcpp)
 sourceCpp('BallMapper.cpp')
 
 //pts <- as.data.frame( read.csv('points.csv',header=F) )
 pts <- rbind( c(0,0), c(1,0), c(10,0), c(12,0), c(30,0), c(33,0) )
 values = as.data.frame( c(1,2,3,4,5,6) )
 //identity orbit, each element is a singleton. 
 orbit = as.data.frame( c(1,2,3,4,5,6) )   
 
 MultiScaleBallMapperGroupActionCpp <- function( points , values , epsilon , orbit )
 {
    output <- SimplifiedMultiScaleBallMapperCppInterfaceGroupActionAndSparseRepresentation( points , values , epsilon , orbit )
    return_list <- output
 }#BallMapperCpp
 
 epsilon = 1
 bm <- MultiScaleBallMapperGroupActionCpp( pts , values , epsilon , orbit )
 BallMapper::ColorIgraphPlot(bm)
 
*/
//This is a version for multiscale ball mapper. The idea is as follows; we select a collection of landmarks on
//some resolution level espilon. We then let the user to manipulate the epsilon to see how the set looks at 
//different resolution levels. 
// [[Rcpp::export]]
List SimplifiedMultiScaleBallMapperCppInterfaceGroupActionAndSparseRepresentation( const DataFrame& points_df , const DataFrame& values_df , double epsilon , const DataFrame& orbit_ )
{
  bool dbg = false                                ;
  if ( points_df.size() == 0 )
  {
    Rcerr << "No points in the SimplifiedMultiScaleBallMapperCppInterfaceGroupActionAndSparseRepresentation procedure, the program will now terminate";
    throw "No points in the SimplifiedMultiScaleBallMapperCppInterfaceGroupActionAndSparseRepresentation procedure, the program will now terminate";
  }

  //the points we obtain from R are unnaturally transposed, so we transpose them back to the situation when we
  //have one point per row.
  //First we need to store them as vector of NumericVectorS:

  std::vector< std::vector< std::pair<unsigned,double> > > points = transpose_points_from_R_to_sparse_vector( points_df );
  
//  Rcerr << "points[0].size()  : " << points[0].size() << endl;

  int number_of_points = points.size();
  if (dbg) Rcerr << "Number of points : " << number_of_points << endl;

  std::vector< std::vector<int> > orbit = transpose_points_from_R_int_version( orbit_ );
  if (dbg) Rcerr << "orbit.size() : " << orbit.size() << endl;


  NumericVector values = values_df[0];
  if (dbg) Rcerr << "values.size() : " << values.size() << endl;


  std::vector< std::vector<size_t> > coverage( number_of_points );
  std::vector< size_t > landmarks;
  landmarks.reserve( (size_t)(0.2*number_of_points) );


  if (dbg) Rcerr << "Entering compute_landmarks_not_transposed_pts_group_action.\n";
  compute_landmarks_not_transposed_pts_group_action_sparse_points( coverage , landmarks ,  points ,  orbit , epsilon );


  if (dbg) Rcerr << "landmarks.size() : " << landmarks.size() << endl;
  //in this case, coverage is not useful, we will need to re-compute it here;
  

  //this is a symetric matrix of a size = number of landmarks. At a (i,j)
  //position of this 
  std::vector< std::vector< double > > graph_incidence_matrix( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )graph_incidence_matrix[i] = std::vector< double >( landmarks.size() );
  
  //now for every landmark:
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
	  //and for every other landmark of a higher index:
	  for ( size_t j = i+1 ; j != landmarks.size() ; ++j )
	  {
		  //for every point in the point cloud, compute its distance to landmarks[i] and landmarks[j], and pick the maximal.
		  //Our purpose is to find a point for which this maximal distance is minimized;
		  double max_distance = std::numeric_limits<double>::max();
		  for ( size_t k = 0 ; k != points.size() ; ++k )
		  {
			  double dist_to_i = compute_distance_standard_points_sparse_points( points[k] , points[ landmarks[i] ] );
			  double dist_to_j = compute_distance_standard_points_sparse_points( points[k] , points[ landmarks[j] ] );
			  double max_distance_from_this_point = std::max( dist_to_i , dist_to_j );
			  if ( max_distance_from_this_point < max_distance )max_distance = max_distance_from_this_point;
		  }
		  graph_incidence_matrix[i][j] = graph_incidence_matrix[j][i] = max_distance;
	  }
  }
  
  //now, for every point we need to points_in_order_from_landmarkscompute how the neigh of points evolve with radius as well as how the coloration is evolving;
  //for every landmark:
  std::vector< std::vector< size_t > > ( landmarks.size() );
  std::vector< std::vector< size_t > > points_in_order_from_landmarks( landmarks.size() );
  std::vector< std::vector< double > > distance_of_points_in_order_from_landmarks( landmarks.size() );
  std::vector< std::vector< double > > coloration_in_order_from_landmarks( landmarks.size() );
  
  
  std::vector< std::vector< std::pair<  double , size_t > > > distance_from_lands_to_points( landmarks.size() );
  //std::vector< std::vector< std::pair< double, double > > > radius_coloration_for_all_landmarks( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
	  std::vector< std::pair<  double , size_t > > distance_point_array( points.size() );
    std::vector< double > distance_point_array_double_component( points.size() );
    std::vector< size_t > distance_point_array_size_t_component( points.size() );
	  //for every point:
	  for ( size_t k = 0 ; k != points.size() ; ++k )
	  {
		  //compute the distance!
		  double dist_to_landmark = compute_distance_standard_points_sparse_points( points[k] , points[ landmarks[i] ] );
		  distance_point_array[k] = std::make_pair( dist_to_landmark , k );
	  }
	  //now we need to sort point_distance_array accordint to 2nd variable:
	  std::sort( distance_point_array.begin() , distance_point_array.end() );
	  for ( size_t aa = 0 ; aa != distance_point_array.size() ; ++aa )
	  {
	    distance_point_array_double_component[aa] = distance_point_array[aa].first;
	    distance_point_array_size_t_component[aa] = distance_point_array[aa].second;
	  }
	  
	  //now we should check how the coloration is evolving as a function of distance from this landmark. 
	  std::vector< double > radius_colororation_for_this_landmark( distance_point_array.size() );
	  double averaged = 0;
	  for ( size_t k = 0 ; k != distance_point_array.size() ; ++k )
	  {
		  averaged *= k;
		  averaged += values[ distance_point_array[ k ].second ];
		  averaged /= (k+1);
		  radius_colororation_for_this_landmark[k] = averaged;
	  }
	  points_in_order_from_landmarks[i] = distance_point_array_size_t_component;
	  distance_of_points_in_order_from_landmarks[i] = distance_point_array_double_component;
	  coloration_in_order_from_landmarks[i] = radius_colororation_for_this_landmark;
  }
  

  NumericVector verts( landmarks.size() );
  for ( size_t i = 0 ; i != landmarks.size() ; ++i )
  {
     verts[i] = i+1;
  }

  List ret;
  ret["vertices"] = verts;
  ret["landmarks"] = landmarks;
  ret["edges"] = graph_incidence_matrix;
  ret["points_in_order_from_landmarks"] = points_in_order_from_landmarks;
  ret["distance_of_points_in_order_from_landmarks"] = distance_of_points_in_order_from_landmarks;
  ret["coloration_in_order_from_landmarks"] = coloration_in_order_from_landmarks;
  return ret;
}



