/*
 *  stdp_connection.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef SGRITTA2017_H
#define SGRITTA2017_H

/* BeginDocumentation
  Name: stdp_synapse - Synapse type for spike-timing dependent
   plasticity.

  Description:
   stdp_synapse is a connector to create synapses with spike time
   dependent plasticity (as defined in [1]). Here the weight dependence
   exponent can be set separately for potentiation and depression.

  Examples:
   multiplicative STDP [2]  mu_plus = mu_minus = 1.0
   additive STDP       [3]  mu_plus = mu_minus = 0.0
   Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
   van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0

  Parameters:
   tau_plus   double - Time constant of STDP window, potentiation in ms
                       (tau_minus defined in post-synaptic neuron)
   lambda     double - Step size
   alpha      double - Asymmetry parameter (scales depressing increments as
                       alpha*lambda)
   mu_plus    double - Weight dependence exponent, potentiation
   mu_minus   double - Weight dependence exponent, depression
   Wmax       double - Maximum allowed weight

  Transmits: SpikeEvent

  References:
   [1] Guetig et al. (2003) Learning Input Correlations through Nonlinear
       Temporally Asymmetric Hebbian Plasticity. Journal of Neuroscience

   [2] Rubin, J., Lee, D. and Sompolinsky, H. (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity, PRL
       86,364-367

   [3] Song, S., Miller, K. D. and Abbott, L. F. (2000). Competitive
       Hebbian learning through spike-timing-dependent synaptic
       plasticity,Nature Neuroscience 3:9,919--926

   [4] van Rossum, M. C. W., Bi, G-Q and Turrigiano, G. G. (2000).
       Stable Hebbian learning from spike timing-dependent
       plasticity, Journal of Neuroscience, 20:23,8812--8821

  FirstVersion: March 2006
  Author: Moritz Helias, Abigail Morrison
  Adapted by: Philipp Weidel
  SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/

// C++ includes:
#include <cmath>
#include <algorithm>
#include <new>
#include <vector>



// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"



std::ofstream amp_;


namespace mynest
{


// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class Sgritta2017 : public nest::Connection< targetidentifierT >
{

public:
  typedef nest::CommonSynapseProperties CommonPropertiesType;
  typedef nest::Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */


  Sgritta2017();

  

  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  Sgritta2017( const Sgritta2017& );

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, nest::ConnectorModel& cm );
 
  
  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  void send( nest::Event& e,
    nest::thread t,
    double t_lastspike,
    const nest::CommonSynapseProperties& cp );

  

  class ConnTestDummyNode : public nest::ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using nest::ConnTestDummyNodeBase::handles_test_event;
    nest::port
    handles_test_event( nest::SpikeEvent&, nest::rport )
    {
      return nest::invalid_port_;
    }
  };

  void
  check_connection( nest::Node& s,
    nest::Node& t,
   nest::rport receptor_type,
    double t_lastspike,
    const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike - get_delay() );
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

private:

 double CalculateMultiplier(double Fpeak){
    double multiplier;
    if(Fpeak < 1 || Fpeak > 10){
      multiplier = 0;
    }

    if(Fpeak <= 6 && Fpeak >= 1){
      multiplier = 0.308*Fpeak - 0.848;
    }
    
    if(Fpeak > 6 && Fpeak <= 10){
      multiplier = -0.115 * Fpeak + 1.69;
    }

    return(multiplier);
  }
     
 double FindPeaks(double perc){

    std::vector<double> DummyFreq;
    std::vector<double> DummyAmp;
    std::vector<double>::iterator found;
    std::vector<double>::iterator max;
    int count = 0;
    double thr;
    double m;
    int index_m;
    double f;
    for(int i = 0; i<W_int;i++){

        if(Frequencies[i] <100){
           DummyAmp.push_back(Amplitudes[i]);
         }
     }

    max = std::max_element(DummyAmp.begin(),DummyAmp.end());
    
    m = *max;
    /*
    thr = m * perc;

    while(flag){
        max = std::max_element(DummyAmp.begin(), DummyAmp.end());
        m = *max;
        found = std::find(DummyAmp.begin(), DummyAmp.end(), m);
        index_m = distance(DummyAmp.begin(), found);

        for(int i = 0; i<= index_m; i++){
        DummyAmp[i] = 0.0;
          }
        f = Frequencies[index_m];
        if(m >= thr){ 
	    DummyFreq.push_back(f);
           }

         if(m < thr){ 
	    flag = 0;
           }

        count++;
    }
     */
    found = std::find(DummyAmp.begin(), DummyAmp.end(), m);
    index_m = distance(DummyAmp.begin(), found);

    return(Frequencies[index_m]);

   }
  

 
 void four1(void)
	{
	    unsigned long n = 0, mmax  = 0, m  = 0, j  = 0, istep  = 0, i  = 0;
	    double wtemp = 0.0, wr = 0.0, wpr = 0.0, wpi = 0.0, wi = 0.0, theta = 0.0;
	    double tempr = 0.0, tempi= 0.0;
	    
	    // reverse-binary reindexing
	    n = W_int<<1;
	    j=1;
            
	    for (i=1; i<n; i+=2) {
		if (j>i) {
		    std::swap(Doppio[j-1], Doppio[i-1]);
		    std::swap(Doppio[j], Doppio[i]);
		}
		m = W_int;
		while (m>=2 && j>m) {
		    j -= m;
		    m >>= 1;
		}
		j += m;
	    };
	 
	    // here begins the Danielson-Lanczos section
	    mmax=2;

	    while (n/2>mmax) {
		istep = mmax<<1;
		theta = -(2*M_PI/mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1; m < mmax; m += 2) {
		    for (i=m; i <= n; i += istep) {
		        j=i+mmax;
		        tempr = wr*Doppio[j-1] - wi*Doppio[j];
		        tempi = wr * Doppio[j] + wi*Doppio[j-1];

		        Doppio[j-1] = Doppio[i-1] - tempr;
		        Doppio[j] = Doppio[i] - tempi;
		        Doppio[i-1] += tempr;
		        Doppio[i] += tempi;
			

		    }
		    wtemp=wr;
		    wr += wr*wpr - wi*wpi;
		    wi += wi*wpr + wtemp*wpi;
		}
		mmax=istep;
	    }



	}

     


   void CalculateA(void){
       double sq_sum;
       double sq_root;
       int j = 0;

       for(int i = 0; i < W_int; i++){
             Amplitudes[i] = 0;
             }

       for(int i=0; i<W_int*2; i = i + 2){
          sq_sum = std::pow(Doppio[i], 2) + std::pow(Doppio[i + 1], 2);
          sq_root = std::sqrt(sq_sum);
          Amplitudes[j] = sq_root;
          j++;
         }

       for(int i = 0; i < W_int; i++){
         if(Frequencies[i] < 1){
             Amplitudes[i] = 0;
             }
       }
    }

   void InstantFreq(double t2, double t1, int P, double A)
      {
      double b = resolution/1000;
      double dT = (t2 - t1)/1000;
      double div = resolution/1000;
      int len = (int)(dT/div + 0.5);
      if (P+1<0 || P+len>=W_int)
                   std::cout << " CHECK7 FAIL " << std::endl;
      for(int i = P + 1; i <= P+len; i++){
         Window[i] = A * std::exp(-1.0*b/0.25);
         b = b + div;
         }

      Window[P+len] = Window[P+len] + 4.0;

       }

   void MoveWindow(double dT, double posOld, int flagM)
    {
     if(flagM == 1 || posOld >= W_int){
        posOld = W_int - 1;}

     int step = dT-((W_int - 1) - posOld);
     if (step<0 || posOld>=W_int)
          std::cout << " CHECK6 FAIL " << " " << step << " " << posOld <<std::endl;
     for(int i = step; i <= posOld; i++)
       {
         Window[i - step] = Window[i];
        }
 
    }

   
   void Inizializza(void){
        double b = 0.0;
        double stepFreq;
        stepFreq = (1000/resolution)/W_int;

        for(int i = 0; i < W_int; i++){  
            Window.push_back(0.0);
            Amplitudes.push_back(0.0);       
            Frequencies.push_back(b);
	    b = b + stepFreq; 

        }
   }   

    
   void Duplica(int flag){

     int j = 0;
     if(flag == 0){
        for(int i = 0; i < W_int*2; i++){
           if(i %2 == 0){
             Doppio.push_back(Window[j]);
             j++;
             
           }
           else if(i %2 != 0){
             Doppio.push_back(0.0);
           }
         
        }

      }
     if(flag != 0){
       for(int i = 0; i < W_int*2; i++){
           if(i %2 == 0){
             Doppio[i] = Window[j];
             j++;
             
           }
           else if(i %2 != 0){
             Doppio[i] = 0.0;
           }
         
        }
     }
   }  

  double
  calculate_k_(double dt)
  {
    double k = 2.0*std::pow(sin(2*M_PI*dt*0.01), 5)*std::exp(-1*std::abs(0.0587701241739 *dt));
   
    return k; 
  }

  double
  facilitate_( double w, double kplus, double scaleFactor, double Peak )
  {
    double norm_w = 0.0;
    double a = 0.0;
    if(Peak >= 1.0 && Peak <= 2.0){
      norm_w = w - std::abs(alpha_*kplus*scaleFactor);
      
     }

     if(Peak > 2.0){
      norm_w = w + (alpha_*kplus*scaleFactor);
      
     }
    return norm_w;
  }

 // double
  //depress_( double w, double kminus )
  //{
    //double norm_w = ( w / Wmax_ )
      //- ( alpha_ * lambda_ * std::pow( w / Wmax_, mu_minus_ ) * kminus );
   // return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
  //}

  // data members of each connection
  double weight_;
  double tau_plus_;
  double lambda_;
  double alpha_;
  long mu_plus_;
  double mu_minus_;
  double Wmax_;
  double Kplus_;
  double Wmin_;
  int p;
  double dtp_;
  double dtn_;
  double t_old;
  double alpha = 0;
  double peak = 0;
  int flag = 0;
  int flagMove = 0;
  int move;
  int posF;
  int pos;
  int pos_old;
  double deltaT;
  int W;
  int W_int;
  double resolution = nest::Time::get_resolution().get_ms();
  std::vector<double> Window; 
  std::vector<double> Frequencies;
  std::vector<double> Doppio;
  std::vector<double> Amplitudes;


  };


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param t_lastspike Time point of last spike emitted
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
Sgritta2017< targetidentifierT >::send( nest::Event& e,
  nest::thread t,
  double t_lastspike,
  const nest::CommonSynapseProperties& )
{
  // synapse STDP depressing/facilitation dynamics
  //   if(t_lastspike >0) {std::cout << "last spike " << t_lastspike <<
  //   std::endl ;}
  double t_spike = e.get_stamp().get_ms();
  // t_lastspike_ = 0 initially

  // use accessor functions (inherited from Connection< >) to obtain delay and
  // target
  nest::Node* target = get_target( t );

  double dendritic_delay = get_delay();

  // get spike history in relevant range (t1, t2] from post-synaptic neuron
  std::deque< nest::histentry >::iterator start;
  std::deque< nest::histentry >::iterator finish;

  // For a new synapse, t_lastspike contains the point in time of the last
  // spike. So we initially read the
  // history(t_last_spike - dendritic_delay, ..., T_spike-dendritic_delay]
  // which increases the access counter for these entries.
  // At registration, all entries' access counters of
  // history[0, ..., t_last_spike - dendritic_delay] have been
  // incremented by Archiving_Node::register_stdp_connection(). See bug #218 for
  // details.
  target->get_history(
    t_lastspike - dendritic_delay, t_spike - dendritic_delay, &start, &finish );
  // facilitation due to post-synaptic spikes since last pre-synaptic spike

 // 

   W = mu_plus_;
   W_int = (int) ((W*1000)/resolution + 0.5);
   
    

    


  while ( start != finish )
  {
    
      posF = (int) (t_spike/resolution + 0.5);
 
      deltaT = (int)((t_spike - t_old)/resolution + 0.5);
      if((pos_old + deltaT) > W_int - 1 && flag != 0 && deltaT < W_int){
          MoveWindow(deltaT, pos_old, flagMove);
          if(flagMove != 0){
               pos = W_int - 1 - deltaT;
	       if (pos<0 || pos>=W_int)
                   std::cout << " CHECK1 FAIL " << std::endl;
               InstantFreq(t_spike, t_old, pos, Window[pos]);           
               Duplica(flagMove); 
               four1();
               CalculateA();
               peak = FindPeaks(mu_minus_);       
               dtp_ = (start ->t_ ) - t_spike; 
	       Kplus_ = calculate_k_(dtp_);
	       alpha = CalculateMultiplier(peak);
	       weight_ = facilitate_( weight_, Kplus_, alpha, peak);
               dtn_ = (start ->t_) - t_lastspike;  
               Kplus_ = calculate_k_(dtn_);
               alpha = CalculateMultiplier(peak);
               weight_ = facilitate_( weight_, Kplus_, alpha, peak);
               if(t == 0){
                   amp_ << FindPeaks(mu_minus_) << "\t" ;}
              }
          if(flagMove == 0){
                p = deltaT - ((W_int  - 1) - pos_old) ; 
                pos = pos_old - p;        
		if (pos<0 || pos>=W_int){
                   std::cout << " CHECK2 FAIL " << pos << std::endl;
                   }
                InstantFreq(t_spike, t_old, pos, Window[pos]);
                Duplica(flagMove);
                four1();
                CalculateA();
                flagMove=1;
                peak = FindPeaks(mu_minus_);
                dtp_ = (start ->t_ ) - t_spike; 
    		Kplus_ = calculate_k_(dtp_);
    		alpha = CalculateMultiplier(peak);
    		weight_ = facilitate_( weight_, Kplus_, alpha, peak);
     	        dtn_ = (start ->t_) - t_lastspike ; 
                Kplus_ = calculate_k_(dtn_);
      	        alpha = CalculateMultiplier(peak);
     	        weight_ = facilitate_( weight_, Kplus_, alpha, peak);
                if(t == 0){
                    amp_ << FindPeaks(mu_minus_) << "\t" ;  }
                }

      }
         
      else if(pos_old + deltaT <= W_int - 1 && flag != 0){
          if (pos_old<0 || pos_old>=W_int)
                   std::cout << " CHECK3 FAIL " << std::endl;  
          InstantFreq(t_spike, t_old, pos_old, Window[pos_old]);
      }
        
      if(flag == 0){
	if(t == 0){
          amp_.close();
          amp_.open( "amp.dat", std::ofstream::app);}
          Inizializza();        
          posF =  posF%W_int;
          if (posF<0 || posF>=W_int)
                   std::cout << " CHECK4 FAIL " << std::endl;
          Window[posF] = 4.0; 
          flag = 1;

      }
      t_old = t_spike;    
      pos_old = posF;
      
	if(t == 0){
        amp_.flush();}


         
       //ciclo if aggiustamento pesi

    start++; // lasciare in fonto qui se no va in vacca tutto il resto;
  }//ciclo while generale
    
	

  // depression due to new pre-synaptic spike
 // weight_ =
   // depress_( weight_, target->get_K_value( t_spike - dendritic_delay ) );

  e.set_receiver( *target );
  if ( weight_ < Wmin_ )
    weight_ = Wmin_;
  if ( weight_ > Wmax_ )
    weight_ = Wmax_;
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
  
  
}


template < typename targetidentifierT >
Sgritta2017< targetidentifierT >::Sgritta2017()
  : ConnectionBase()
  , weight_( 1.0 )
  , tau_plus_( 20.0 )
  , lambda_( 0.01 )
  , alpha_( 1.0 )
  , mu_plus_( 1 )
  , mu_minus_( 1.0 )
  , Wmax_( 100.0 )
  , Kplus_( 0.0 )
  , Wmin_(-100.0)
{
}

template < typename targetidentifierT >
Sgritta2017< targetidentifierT >::Sgritta2017(
  const Sgritta2017< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , tau_plus_( rhs.tau_plus_ )
  , lambda_( rhs.lambda_ )
  , alpha_( rhs.alpha_ )
  , mu_plus_( rhs.mu_plus_ )
  , mu_minus_( rhs.mu_minus_ )
  , Wmax_( rhs.Wmax_ )
  , Wmin_( rhs.Wmin_ )
  , Kplus_( rhs.Kplus_ )
{
}

template < typename targetidentifierT >
void
Sgritta2017< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, nest::names::weight, weight_ );
  def< double >( d, nest::names::tau_plus, tau_plus_ );
  def< double >( d, nest::names::lambda, lambda_ );
  def< double >( d, nest::names::alpha, alpha_ );
  def< long >( d, nest::names::mu_plus, mu_plus_ );
  def< double >( d, nest::names::mu_minus, mu_minus_ );
  def< double >( d, nest::names::Wmax, Wmax_ );
  def< double >( d, nest::names::Wmin, Wmin_ );
  def< long >( d, nest::names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
Sgritta2017< targetidentifierT >::set_status( const DictionaryDatum& d,
  nest::ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, nest::names::weight, weight_ );
  updateValue< double >( d, nest::names::tau_plus, tau_plus_ );
  updateValue< double >( d, nest::names::lambda, lambda_ );
  updateValue< double >( d, nest::names::alpha, alpha_ );
  updateValue< long  >( d, nest::names::mu_plus, mu_plus_ );
  updateValue< double >( d, nest::names::mu_minus, mu_minus_ );
  updateValue< double >( d, nest::names::Wmax, Wmax_ );
  updateValue< double >( d, nest::names::Wmin, Wmin_ );
  amp_ <<std::endl;
  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) )
         == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw nest::BadProperty( "Weight aed Wmax must have same sign." );
  }
}

} // of namespace nest

#endif // of #ifndef SGRITTA2017_Hhttps://open.spotify.com/track/1mBMgkXH9n2SUQSB4LmihW