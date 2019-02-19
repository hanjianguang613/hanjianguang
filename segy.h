#include "gb.h"

#ifndef SEGY_H_
#define SEGY_H_

#define SEGY_ASCII_HDR_SIZE 3200     //Ascii Volumn Header 
#define SEGY_VOLUMN_HDR_SIZE 400     //Binary Volumn Header
#define SEGY_TRACE_HDR_SIZE 240      //Trace Header

#ifndef UINT8_T
#define UINT8_T
typedef unsigned char uint8_t;
#endif

#ifndef INT16_T
#define INT16_T
typedef short int16_t;
#endif

#ifndef INT32_T
#define INT32_T
typedef long int32_T;
#endif


//--------Binary Volumn Header-------//
typedef struct segy_volumn_hdrstruct
{
	int32_T job_id_number;						 // 3201-3204   4    
	int32_T line_number;						 // 3205-3208   4   
	int32_T reel_number;						 // 3209-3212   4    
	int16_t traces_per_record;					 // 3213-3214   2=   
	int16_t aux_traces_per_record;			  	 // 3215-3216   2    
	int16_t sample_data_interval_ms;			 // 3217-3218   2=   
	int16_t original_data_interval_ms;			 // 3219-3220   2    
	int16_t samples_per_trace;					 // 3221-3222   2=   
	int16_t original_samples_per_trace;			 // 3223-3224   2    
	int16_t data_sample_format_code;			 // 3225-3226   2=   
	int16_t CDP_fold;							 // 3227-3228   2=   
	int16_t trace_sorting_code;					 // 3229-3230   2=   
	int16_t vertical_sum_code;					 // 3231-3232   2    
	int16_t sweep_frequency_start_hz;			 // 3233-3234   2    
	int16_t sweep_frequency_end_hz;				 // 3235-3236   2    
	int16_t sweep_length_ms;					 // 3237-3238   2   
	int16_t sweep_type_code;					 // 3239-3240   2    
	int16_t trace_number_of_sweep_channel;       // 3241-3242   2    
	int16_t sweep_trace_taper_length_start_ms;   // 3243-3244   2    
	int16_t sweep_trace_taper_length_end_ms;     // 3245-3246   2    
	int16_t taper_type_code;                     // 3247-3248   2    
	int16_t correlated_data_traces_flag;         // 3249-3250   2    
	int16_t binary_gain_recovered_flag;          // 3251-3252   2    
	int16_t amplitude_recovery_method_code;      // 3253-3254   2    
	int16_t measurement_system;                  // 3255-3256   2=    
	int16_t impulse_signal_polarity;             // 3257-3258   2    
	int16_t vibratory_polarity_code;             // 3259-3260   2    
	uint8_t buffer[340];						 // 12 bytes + 12 bytes + 36 bytes + 340 buffer bytes = REEL_HDR_SIZE  // 3261-3600   340  

} Segy_volumn_hdr;

//--------Trace Header----------//
typedef struct segy_trace_hdrstruct 
{
	int32_T trace_sequence_number_within_line;							 // 1-4     4   
	int32_T trace_sequence_number_within_reel;							 // 5-8     4  
	int32_T original_field_record_number;								 // 9-12    4*  
	int32_T trace_sequence_number_within_original_field_record;			 // 13-16   4  
	int32_T energy_source_point_number;									 // 17-20   4   
	int32_T cdp_ensemble_number;										 // 21-24   4* 
	int32_T trace_sequence_number_within_cdp_ensemble;					 // 25-28   4*  
	int16_t trace_identification_code;									 // 29-30   2   
	int16_t number_of_vertically_summed_traces_yielding_this_trace;      // 31-32   2   
	int16_t number_of_horizontally_stacked_traced_yielding_this_trace;   // 33-34   2   
	int16_t data_use;													 // 35-36   2   
	int32_T distance_from_source_point_to_receiver_group;				 // 37-40   4*  
	int32_T receiver_group_elevation;									 // 41-44   4   
	int32_T surface_elevation_at_source;								 // 45-48   4   
	int32_T source_depth_below_surface;									 // 49-52   4   
	int32_T datum_elevation_at_receiver_group;							 // 53-56   4   
	int32_T datum_elevation_at_source;									 // 57-60   4   
	int32_T water_depth_at_source;										 // 61-64   4   
	int32_T water_depth_at_receiver_group;								 // 65-68   4   
	int16_t scalar_for_elevations_and_depths;							 // 69-70   2   
	int16_t scalar_for_coordinates;										 // 71-72   2   
	int32_T x_source_coordinate;										 // 73-76   4   
	int32_T y_source_coordinate;										 // 77-80   4   
	int32_T x_receiver_group_coordinate;								 // 81-84   4  
	int32_T y_receiver_group_coordinate;							     // 85-88   4   
	int16_t coordinate_units;											 // 89-90   2   

	int16_t weathering_velocity;										 // 91-92   2   
	int16_t subweathering_velocity;										 // 93-94   2   
	int16_t uphole_time_at_source;									     // 95-96   2   
	int16_t uphole_time_at_group;										 // 97-98   2   
	int16_t source_static_correction;									 // 99-100  2   

	int16_t group_static_correction;            // 101-102 2   
	int16_t total_static_applied;               // 103-104 2   
	int16_t lag_time_a;                         // 105-106 2   
	int16_t lag_time_b;                         // 107-108 2  
	int16_t delay_according_time;               // 109-110 2   

	int16_t mute_time_start;                    // 111-112 2   
	int16_t mute_time_end;                      // 113-114 2   
	int16_t samples_in_this_trace;              // 115-116 2* 
	int16_t sample_intervall;                   // 117-118 2* 
	int16_t gain_type_instruments;              // 119-120 2   

	int16_t igc;	           // 121-122 2   
	int16_t igi;		       // 123-124 2   
	int16_t corr;	           // 125-126 2   
	int16_t sfs;	           // 127-128 2   
	int16_t sfe;	           // 129-130 2   
	int16_t slen;              // 131-132 2   
	int16_t styp;	           // 133-134 2   
	int16_t stas;	           // 135-136 2   
	int16_t stae;	           // 137-138 2   
	int16_t tatyp;             // 139-140 2   
	int16_t afilf;			   // 141-142 2   
	int16_t afils;	           // 143-144 2   
	int16_t nofilf;	           // 145-146 2   
	int16_t nofils;	           // 147-148 2   
	int16_t lcf;	           // 149-150 2   
	int16_t hcf;	           // 151-152 2   
	int16_t lcs;	           // 153-154 2   
	int16_t hcs;	           // 155-156 2   
	int16_t year;	           // 157-158 2  
	int16_t day;	           // 159-160 2   
	int16_t hour;	           // 161-162 2   
	int16_t minute;	           // 163-164 2   
	int16_t sec;	           // 165-166 2   
	int16_t timbas;	           // 167-168 2	  I = local.	2 = GMT. 3 = other.
	int16_t trwf;	           // 169-170 2   trace weighting factor
	int16_t grnors;	           // 171-172 2   geophone group number of roll switch position one
	int16_t grnofr;	           // 173-174 2   geophone group number of trace one within original field record
	int16_t grnlof;	           // 175-176 2   geophone group number of last trace within original field record
	int16_t gaps;	           // 177-178 2   gap size (total number of groups dropped)
	int16_t otrav;	           // 179-180 2   overtravel taper code: 1 = down (or behind) 2 = up (or ahead)
	int32_T user_define[15];   // 181-240 60  

} Segy_trace_hdr;

	/* read the volumn header information */
	Segy_volumn_hdr Readvolumn(char *infile);
    /*Calculate how many traces in this record*/
	int get_trace_num(int samples_per_trace, char *infile);
	/* read segy file */
	void Readsegy(int trace_num, int samples_per_trace, Segy_trace_hdr *trace_hdr,float **trace_data, char *infile);
	/* write segy file */
	void Writesegy(int trace_num, int samples_per_trace,Segy_volumn_hdr volumn_hdr, Segy_trace_hdr *trace_hdr,float **trace_data, char *outfile);
	
#endif