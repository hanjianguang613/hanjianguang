/*****************************************************
Structure:
typedef struct segy_volumn_hdr 
{
	int32_T job_id_number;              
	int32_T line_number;                
	int32_T reel_number;                
	int16_t traces_per_record;         
	int16_t aux_traces_per_record;      
	int16_t sample_data_interval_ms;   
	int16_t original_data_interval_ms; 
	int16_t samples_per_trace;         
	int16_t original_samples_per_trace; 
	int16_t data_sample_format_code;    
	int16_t CDP_fold;                   
	int16_t trace_sorting_code;        
	int16_t vertical_sum_code;          
	int16_t sweep_frequency_start_hz;  
	int16_t sweep_frequency_end_hz;    
	int16_t sweep_length_ms;          
	int16_t sweep_type_code;            
	int16_t trace_number_of_sweep_channel;      
	int16_t sweep_trace_taper_length_start_ms; 
	int16_t sweep_trace_taper_length_end_ms;    
	int16_t taper_type_code;                   
	int16_t correlated_data_traces_flag;       
	int16_t binary_gain_recovered_flag;        
	int16_t amplitude_recovery_method_code;     
	int16_t measurement_system;                 
	int16_t impulse_signal_polarity;            
	int16_t vibratory_polarity_code;           
	uint8_t buffer[340]; 
} segy_volumn_hdr;

//--------Trace Header----------//
typedef struct segy_trace_hdr 
{
	int32_T trace_sequence_number_within_line;  
	int32_T trace_sequence_number_within_reel; 
	int32_T original_field_record_number;       
	int32_T trace_sequence_number_within_original_field_record; 
	int32_T energy_source_point_number;        
	int32_T cdp_ensemble_number;              
	int32_T trace_sequence_number_within_cdp_ensemble;     
	int16_t trace_identification_code;        
	int16_t number_of_vertically_summed_traces_yielding_this_trace;   
	int16_t number_of_horizontally_stacked_traced_yielding_this_trace;  
	int16_t data_use;                         
	int32_T distance_from_source_point_to_receiver_group;   
	int32_T receiver_group_elevation;        
	int32_T surface_elevation_at_source;      
	int32_T source_depth_below_surface;     
	int32_T datum_elevation_at_receiver_group; 
	int32_T datum_elevation_at_source;        
	int32_T water_depth_at_source;           
	int32_T water_depth_at_receiver_group;  
	int16_t scalar_for_elevations_and_depths;   
	int16_t scalar_for_coordinates;            
	int32_T x_source_coordinate;               
	int32_T y_source_coordinate;               
	int32_T x_receiver_group_coordinate;       
	int32_T y_receiver_group_coordinate;       
	int16_t coordinate_units;                 

	int16_t weathering_velocity;               
	int16_t subweathering_velocity;           
	int16_t uphole_time_at_source;             
	int16_t uphole_time_at_group;            
	int16_t source_static_correction;        

	int16_t group_static_correction;           
	int16_t total_static_applied;              
	int16_t lag_time_a;                        
	int16_t lag_time_b;                         
	int16_t delay_according_time;           

	int16_t mute_time_start;               
	int16_t mute_time_end;         
	int16_t samples_in_this_trace;         
	int16_t sample_intervall;       
	int16_t gain_type_instruments;      

	int16_t igc;	
	int16_t igi;		
	int16_t corr;	    
	int16_t sfs;	   
	int16_t sfe;	  
	int16_t slen;       
	int16_t styp;	  
	int16_t stas;	      
	int16_t stae;	  
	int16_t tatyp;  
	int16_t afilf;	   
	int16_t afils;	   
	int16_t nofilf;	    
	int16_t nofils;	   
	int16_t lcf;	    
	int16_t hcf;	        
	int16_t lcs;	      
	int16_t hcs;	   
	int16_t year;	      
	int16_t day;	
	int16_t hour;	   
	int16_t minute;	    
	int16_t sec;	       
	int16_t timbas;	  
	int16_t trwf;	      
	int16_t grnors;	  
	int16_t grnofr;	     
	int16_t grnlof;	    
	int16_t gaps;	     
	int16_t otrav;	     
	int32_T user_define[15];

} segy_trace_hdr;

*********************************************************************************************
Function Prototypes:
	segy_volumn_hdr *Readvolumn(char *infile);
	int get_trace_num(int samples_per_trace, char *infile);
	void Readsegy(int trace_num, int samples_per_trace, segy_trace_hdr *trace_hdr,float **trace_data, char *infile);
	void Writesegy(int trace_num, int samples_per_trace,segy_volumn_hdr *volumn_hdr, segy_trace_hdr *trace_hdr,float **trace_data, char *outfile);
	
********************************************************************************************/

/* ******************************************************************************************
 Author: Han Jianguang
 Note: The program achieved read and write segy file.
*********************************************************************************************/

#include "segy.h"

/* read the volumn header information */
Segy_volumn_hdr Readvolumn(char *infile)
{
	Segy_volumn_hdr volumn_hdr;
	int i;
	FILE *fp1=NULL;
	fp1=fopen(infile,"rb");
	if(fp1==NULL)
	{
		printf("hjg.segy is not exist !\n");
		return ;
	}
	fseek(fp1,3200L,0);
	/* Read Binary Volumn Header */
	fread(&volumn_hdr.job_id_number, sizeof(int32_T),1,fp1);
	fread(&volumn_hdr.line_number, sizeof(int32_T),1,fp1);
	fread(&volumn_hdr.reel_number, sizeof(int32_T),1,fp1);
	fread(&volumn_hdr.traces_per_record, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.aux_traces_per_record, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.sample_data_interval_ms, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.original_data_interval_ms, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.samples_per_trace, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.original_samples_per_trace, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.data_sample_format_code, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.CDP_fold, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.trace_sorting_code, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.vertical_sum_code, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.sweep_frequency_start_hz, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.sweep_frequency_end_hz, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.sweep_length_ms, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.sweep_type_code, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.trace_number_of_sweep_channel, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.sweep_trace_taper_length_start_ms, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.sweep_trace_taper_length_end_ms, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.taper_type_code, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.correlated_data_traces_flag, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.binary_gain_recovered_flag, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.amplitude_recovery_method_code, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.measurement_system, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.impulse_signal_polarity, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.vibratory_polarity_code, sizeof(int16_t),1,fp1);
	fread(&volumn_hdr.buffer, sizeof(uint8_t),340,fp1);

	fclose(fp1);
	return volumn_hdr;
	
}

/*Calculate how many traces in this record*/
int get_trace_num(int samples_per_trace, char *infile)
{
	int trace_num;
	long pos;
	FILE *fp2=NULL;
	fp2=fopen(infile,"rb");
	fseek(fp2,0L,2);
	pos=ftell(fp2);

	trace_num = (pos - SEGY_ASCII_HDR_SIZE - SEGY_VOLUMN_HDR_SIZE) / (samples_per_trace *4 + SEGY_TRACE_HDR_SIZE);
	fclose(fp2);
	return trace_num;
}

/* Read Trace Header and Trace Data */
void Readsegy(int trace_num, int samples_per_trace, Segy_trace_hdr *trace_hdr,float **trace_data, char *infile)
{
	int i;
	FILE *fp3=NULL;
	fp3=fopen(infile,"rb");
	fseek(fp3,3600L,0);
	/* Read Trace Header */
	for(i=0;i<trace_num;i++)
	{
		fread(&trace_hdr[i].trace_sequence_number_within_line, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].trace_sequence_number_within_reel, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].original_field_record_number, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].trace_sequence_number_within_original_field_record, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].energy_source_point_number, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].cdp_ensemble_number, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].trace_sequence_number_within_cdp_ensemble, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].trace_identification_code, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].number_of_vertically_summed_traces_yielding_this_trace, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].number_of_horizontally_stacked_traced_yielding_this_trace, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].data_use, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].distance_from_source_point_to_receiver_group, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].receiver_group_elevation, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].surface_elevation_at_source, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].source_depth_below_surface, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].datum_elevation_at_receiver_group, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].datum_elevation_at_source, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].water_depth_at_source, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].water_depth_at_receiver_group, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].scalar_for_elevations_and_depths, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].scalar_for_coordinates, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].x_source_coordinate, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].y_source_coordinate, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].x_receiver_group_coordinate, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].y_receiver_group_coordinate, sizeof(int32_T),1,fp3);
		fread(&trace_hdr[i].coordinate_units, sizeof(int16_t),1,fp3);

		fread(&trace_hdr[i].weathering_velocity, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].subweathering_velocity, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].uphole_time_at_source, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].uphole_time_at_group, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].source_static_correction, sizeof(int16_t),1,fp3);

		fread(&trace_hdr[i].group_static_correction, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].total_static_applied, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].lag_time_a, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].lag_time_b, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].delay_according_time, sizeof(int16_t),1,fp3);

		fread(&trace_hdr[i].mute_time_start, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].mute_time_end, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].samples_in_this_trace, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].sample_intervall, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].gain_type_instruments, sizeof(int16_t),1,fp3);

		fread(&trace_hdr[i].igc, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].igi, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].corr, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].sfs, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].sfe, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].slen, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].styp, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].stas, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].stae, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].tatyp, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].afilf, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].afils, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].nofilf, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].nofils, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].lcf, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].hcf, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].lcs, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].hcs, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].year, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].day, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].hour, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].minute, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].sec, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].timbas, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].trwf, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].grnors, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].grnofr, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].grnlof, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].gaps, sizeof(int16_t),1,fp3);
		fread(&trace_hdr[i].otrav, sizeof(int16_t),1,fp3);

		fread(trace_hdr[i].user_define, sizeof(int32_T),15,fp3);
		/* Read Trace data */
		fread(trace_data[i],sizeof(float),samples_per_trace,fp3);
	}
	fclose(fp3);
}

/* write segy file */
void Writesegy(int trace_num, int samples_per_trace,Segy_volumn_hdr volumn_hdr, Segy_trace_hdr *trace_hdr,float **trace_data, char *outfile)
{
	char wjt[SEGY_ASCII_HDR_SIZE];
	int i;
	FILE *fp4=NULL;
	fp4=fopen(outfile,"wb");
	/* Write Ascii Header */
	fwrite(wjt,1,3200,fp4);
	/* Read Binary Volumn Header */
	fwrite(&volumn_hdr.job_id_number, sizeof(int32_T),1,fp4);
	fwrite(&volumn_hdr.line_number, sizeof(int32_T),1,fp4);
	fwrite(&volumn_hdr.reel_number, sizeof(int32_T),1,fp4);
	fwrite(&volumn_hdr.traces_per_record, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.aux_traces_per_record, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.sample_data_interval_ms, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.original_data_interval_ms, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.samples_per_trace, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.original_samples_per_trace, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.data_sample_format_code, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.CDP_fold, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.trace_sorting_code, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.vertical_sum_code, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.sweep_frequency_start_hz, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.sweep_frequency_end_hz, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.sweep_length_ms, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.sweep_type_code, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.trace_number_of_sweep_channel, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.sweep_trace_taper_length_start_ms, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.sweep_trace_taper_length_end_ms, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.taper_type_code, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.correlated_data_traces_flag, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.binary_gain_recovered_flag, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.amplitude_recovery_method_code, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.measurement_system, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.impulse_signal_polarity, sizeof(int16_t),1,fp4);
	fwrite(&volumn_hdr.vibratory_polarity_code, sizeof(int16_t),1,fp4);
	fwrite(volumn_hdr.buffer, sizeof(uint8_t),340,fp4);

	/*Write Trace Header */
	for(i=0;i<trace_num;i++)
	{
		fwrite(&trace_hdr[i].trace_sequence_number_within_line, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].trace_sequence_number_within_reel, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].original_field_record_number, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].trace_sequence_number_within_original_field_record, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].energy_source_point_number, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].cdp_ensemble_number, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].trace_sequence_number_within_cdp_ensemble, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].trace_identification_code, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].number_of_vertically_summed_traces_yielding_this_trace, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].number_of_horizontally_stacked_traced_yielding_this_trace, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].data_use, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].distance_from_source_point_to_receiver_group, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].receiver_group_elevation, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].surface_elevation_at_source, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].source_depth_below_surface, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].datum_elevation_at_receiver_group, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].datum_elevation_at_source, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].water_depth_at_source, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].water_depth_at_receiver_group, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].scalar_for_elevations_and_depths, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].scalar_for_coordinates, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].x_source_coordinate, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].y_source_coordinate, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].x_receiver_group_coordinate, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].y_receiver_group_coordinate, sizeof(int32_T),1,fp4);
		fwrite(&trace_hdr[i].coordinate_units, sizeof(int16_t),1,fp4);

		fwrite(&trace_hdr[i].weathering_velocity, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].subweathering_velocity, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].uphole_time_at_source, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].uphole_time_at_group, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].source_static_correction, sizeof(int16_t),1,fp4);

		fwrite(&trace_hdr[i].group_static_correction, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].total_static_applied, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].lag_time_a, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].lag_time_b, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].delay_according_time, sizeof(int16_t),1,fp4);

		fwrite(&trace_hdr[i].mute_time_start, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].mute_time_end, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].samples_in_this_trace, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].sample_intervall, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].gain_type_instruments, sizeof(int16_t),1,fp4);

		fwrite(&trace_hdr[i].igc, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].igi, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].corr, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].sfs, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].sfe, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].slen, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].styp, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].stas, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].stae, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].tatyp, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].afilf, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].afils, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].nofilf, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].nofils, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].lcf, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].hcf, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].lcs, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].hcs, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].year, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].day, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].hour, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].minute, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].sec, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].timbas, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].trwf, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].grnors, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].grnofr, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].grnlof, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].gaps, sizeof(int16_t),1,fp4);
		fwrite(&trace_hdr[i].otrav, sizeof(int16_t),1,fp4);

		fwrite(trace_hdr[i].user_define, sizeof(int32_T),15,fp4);
		/* Write Trace data */
		fwrite(trace_data[i],sizeof(float),samples_per_trace,fp4);
	}

	fclose(fp4);

}
