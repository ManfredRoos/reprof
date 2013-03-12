#include "includes/fann.h"
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include "util.h"
#include "reprof_struct.h"
#include "blastpsimat.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//
//  main.c
//  madprof
// #--------------------------------------------------
// Predict secondary structure and solvent
// accessibility from sequece
//
// steinegger_martin@web.de
//  Created by Martin Steinegger on 25.03.12.
//  Copyright (c) 2012 -. All rights reserved.
//
// [MR] added some bugfixes, mainly to make the Code more readable, function wrappers to avoid typos, some less copy overhead.

/**function declarations*/
features_list * parse_feature_file(const char *configuration_file);
int parse_args (int argc, const char * argv[]);

void help();
int build_models_and_features();
void set_snap();
int f_exists();
void free_oris();
void free_models_and_features();
void Precompute();
void do_the_expensive_stuff();
const char *get_filename_ext(const char *filepath);

 char *get_filename( char *filepath);
float_array_2d * create_and_run(int Model_num ,float_array_2d * pre_output);
float_array_2d * create_and_run_from_to(int Model_num ,float_array_2d * pre_output,size_t from,size_t to);
float_array_2d * copy(float_array_2d * _ori);

/**Global variable definitions*/
char * pwd;
//int exit=-1;
char * in_file=NULL;
char * out_file=NULL;
char * sequ=NULL;
const char * pathToModel="/var/tmp/mamut/Mprof/"; //→not in userspace   "/usr/share/madprof/";
const char * pathToOut=NULL;							// «mr   // use the NULL not 0
int hasPsiMat=0;
char outPahtBuf[256];
int max_window=0;
size_t input_size=0;
size_t chain_length=0;


int snap_20=20;				// «mr   // make 19 non native:0 or 20 including  the native :1
int dbug=0;
int short_output=0;
int quiet=0;

struct fann **models=NULL;
features_list **features=NULL;

enum
   {
	FA_MODEL 	= 0,
	FU_MODEL 	= 1,
	FB_MODEL 	= 2,
	FUU_MODEL 	= 3,
	FUB_MODEL 	= 4,
	FBU_MODEL 	= 5,
	FBB_MODEL 	= 6,
	A_MODEL 	= 7,
	U_MODEL 	= 8,
	B_MODEL 	= 9,
	UU_MODEL 	= 10,
	UB_MODEL 	= 11,
	BU_MODEL 	= 12,
	BB_MODEL 	= 13,
	model_size 	= 14
   } Models;
//const int FA_MODEL = 0;
//const int FU_MODEL = 1;
//const int FB_MODEL = 2;
//const int FUU_MODEL = 3;
//const int FUB_MODEL = 4;
//const int FBU_MODEL = 5;
//const int FBB_MODEL = 6;
//const int A_MODEL = 7;
//const int U_MODEL = 8;
//const int B_MODEL = 9;
//const int UU_MODEL = 10;
//const int UB_MODEL = 11;
//const int BU_MODEL = 12;
//const int BB_MODEL = 13;
//   const int model_size = 14;

const char model_names[14][4] = {"fa","fu", "fb", "fuu", "fub", "fbu", "fbb", "a", "u", "b", "uu", "ub", "bu", "bb"};

float_array_2d * a, * u, * b,  * uu,* ub ,* bu,* bb,* sec ;
float_array_2d * sec_ori=NULL,* a_ori=NULL,* u_ori=NULL,* b_ori=NULL;
float_array_2d * sec_tmp, * a_tmp,* b_tmp ,* u_tmp;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** short main to make readable whats happening */
int main (int argc, const char * argv[])
{	
  
    pwd=get_current_dir_name(); //  getcwd(pwd,sizeof(pwd));		// where are we?

    parse_args(argc, argv); /// faild to find -i  → -13 , -o → -14
       
    build_models_and_features(); /// failed to build models → -12
    
    Precompute(); 

    do_the_expensive_stuff();
    
    ///////// free up before the end
    
    free_oris();
        
    free_models_and_features();

    return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/** Load models and feature lists*/
int build_models_and_features(){
  if(!f_exists(pathToModel ) ){
    printf("model dir [%s] not found ... maybe try \" --model /var/tmp/mamut/Mprof_models/ \"   \n", pathToModel);
    exit( -13);// return -12;
  }
    models = malloc(model_size * sizeof(struct fann* ));

    features = malloc(model_size * sizeof(struct features_in_out* ));

    for(int i  = 0; i<model_size; i++){
        char modelPahtBuf[256];
        char featurePahtBuf[256];
        snprintf(modelPahtBuf, sizeof modelPahtBuf, "%s%s%s%s", pathToModel,"//",(char*) &model_names[i],".model");
        snprintf(featurePahtBuf, sizeof featurePahtBuf, "%s%s%s%s", pathToModel,"//",(char*) &model_names[i],".features");
        models[i] = fann_create_from_file(modelPahtBuf);
        features[i] = parse_feature_file(featurePahtBuf);
    }
    if(quiet<1)  printf("models read\n");
     return 0;
}

void free_oris(){
    free_float_array_2d(a_ori);							//« mr : was missing so we allocated a lot of ram over time without ever freeing it
    free_float_array_2d(u_ori);
    free_float_array_2d(b_ori);
    free_float_array_2d(sec_ori);
  
}


void free_models_and_features(){

for(int i = 0; i < sizeof(models) / sizeof(struct fann); i++){
    free(models[i]);
}
for(int i = 0; i < sizeof(features) / sizeof(features_list); i++){
    free_feature_list(features[i]);
}
}
/**go once trough the input and predict everythign for the default sequence*/
void Precompute(){
	 //Precompute
	    //Do prediction for the inputfile
	 (quiet<1||dbug>0) &&    printf("seq main: %s %lu\n", sequ,chain_length);
	    if(hasPsiMat==1){
	 quiet<1 &&   printf("hasPsiMat\n");										//»mr dbug
	         input_size=(chain_length - 1) - 0+1;
	         a=create_and_run(A_MODEL,NULL);
	         b=create_and_run(B_MODEL,NULL);
	         u=create_and_run(U_MODEL,NULL);

	//         a = run_model(models[A_MODEL], create_inputs(features[A_MODEL], 0, chain_length - 1,sequ,chain_length,models[A_MODEL]->num_input,NULL),input_size);
	//         u = run_model(models[U_MODEL], create_inputs(features[U_MODEL], 0, chain_length - 1,sequ,chain_length,models[U_MODEL]->num_input,NULL),input_size);
	//         b = run_model(models[B_MODEL], create_inputs(features[B_MODEL], 0, chain_length - 1,sequ,chain_length,models[B_MODEL]->num_input,NULL),input_size);
	//
	         uu=create_and_run(UU_MODEL,u);
	         ub=create_and_run(UB_MODEL,u);
	         bu=create_and_run(BU_MODEL,b);
	         bb=create_and_run(BB_MODEL,b);
	//         uu = run_model(models[UU_MODEL], create_inputs(features[UU_MODEL], 0, chain_length - 1,sequ,chain_length,models[UU_MODEL]->num_input, u),input_size);
	//         ub = run_model(models[UB_MODEL], create_inputs(features[UB_MODEL], 0, chain_length - 1,sequ,chain_length,models[UB_MODEL]->num_input, u),input_size);
	//         bu = run_model(models[BU_MODEL], create_inputs(features[BU_MODEL], 0, chain_length - 1,sequ,chain_length,models[BU_MODEL]->num_input, b),input_size);
	//         bb = run_model(models[BB_MODEL], create_inputs(features[BB_MODEL], 0, chain_length - 1,sequ,chain_length,models[BB_MODEL]->num_input, b),input_size);
	//
	         sec_ori=jury(uu, ub, bu, bb);

	        write_output(&outPahtBuf, sec_ori, a, sequ, chain_length);  					//← code duplication result  will be lost by overwrite →↓
	        free_float_array_2d(a);
	        free_float_array_2d(u);
	        free_float_array_2d(b);
	        free_float_array_2d(uu);
	        free_float_array_2d(ub);
	        free_float_array_2d(bu);
	        free_float_array_2d(bb);
	        free_float_array_2d(sec_ori);
	    }
	   (quiet<2||dbug>0) &&    printf("input_size: %d\n",chain_length);									//»mr dbug
	     input_size=(chain_length - 1) - 0+1;
	   if(quiet<1)       printf("run:\n");														//»mr dbug
	        a_ori=create_and_run(FA_MODEL,NULL);
	        u_ori=create_and_run(FU_MODEL,NULL);
	        b_ori=create_and_run(FB_MODEL,NULL);
	//    float_array_2d * a_ori = run_model(models[FA_MODEL], create_inputs(features[FA_MODEL], 0, chain_length - 1,sequ,chain_length,models[FA_MODEL]->num_input,NULL),input_size);
	//    float_array_2d * u_ori = run_model(models[FU_MODEL], create_inputs(features[FU_MODEL], 0, chain_length - 1,sequ,chain_length,models[FU_MODEL]->num_input,NULL),input_size);
	//    float_array_2d * b_ori = run_model(models[FB_MODEL], create_inputs(features[FB_MODEL], 0, chain_length - 1,sequ,chain_length,models[FB_MODEL]->num_input,NULL),input_size);
	   if(quiet<1)       printf("run:\n");									//»mr dbug
	        uu=create_and_run(FUU_MODEL,u_ori);
	        ub=create_and_run(FUB_MODEL,u_ori);
	        bu=create_and_run(FBU_MODEL,b_ori);
	        bb=create_and_run(FBB_MODEL,b_ori);
	//    float_array_2d * uu = run_model(models[FUU_MODEL], create_inputs(features[FUU_MODEL], 0, chain_length - 1,sequ,chain_length,models[FUU_MODEL]->num_input, u_ori),input_size);
	//    float_array_2d * ub = run_model(models[FUB_MODEL], create_inputs(features[FUB_MODEL], 0, chain_length - 1,sequ,chain_length,models[FUB_MODEL]->num_input, u_ori),input_size);
	//    float_array_2d * bu = run_model(models[FBU_MODEL], create_inputs(features[FBU_MODEL], 0, chain_length - 1,sequ,chain_length,models[FBU_MODEL]->num_input, b_ori),input_size);
	//    float_array_2d * bb = run_model(models[FBB_MODEL], create_inputs(features[FBB_MODEL], 0, chain_length - 1,sequ,chain_length,models[FBB_MODEL]->num_input, b_ori),input_size);
	//
	//    float_array_2d *
	    sec_ori=jury(uu, ub, bu, bb);
	    char oriPahtBuf[256];
	  if(quiet<1)     printf("run:\n");															//»mr dbug
	    snprintf(oriPahtBuf, sizeof oriPahtBuf, "%s_ORI", outPahtBuf);
	    write_output(&oriPahtBuf, sec_ori, a_ori, sequ, chain_length);                       //← code duplication result(file) of above will be overwritten!!, maybe else{} missing?
	    //Free memory
	    free_float_array_2d(uu);
	    free_float_array_2d(ub);
	    free_float_array_2d(bu);
	    free_float_array_2d(bb);


}

/** the heavy duty computation is in here: for all mutations predict the +- maxwindow arround it and it selve, and integrate that in the original prediction to write a consistant outpuit file*/
void do_the_expensive_stuff(){
  if(quiet<1)    printf("run2: the really heavy stuff:\n");								//»mr dbug
    if(sec_ori==NULL||a_ori==NULL||u_ori==NULL||b_ori==NULL)exit( 3);
    char extra_out[256];
    snprintf(extra_out, sizeof extra_out, "%s%s", pathToOut,"_EXTRA_");
    FILE *EXTRA = fopen(extra_out, "w");
    
   for(int pos = 0; pos < chain_length; pos++){
       char aa_ori = sequ[pos];
       size_t seq_from = max(((pos - 1) - 1 * max_window + 1), 0);
       size_t seq_to =   min(((pos - 1) + 1 * max_window - 1), (chain_length - 1));

       size_t struc_from = max(((pos - 1) - 2 * max_window + 2), 0);
       size_t struc_to   = min(((pos - 1) + 2 * max_window - 2), (chain_length - 1));

/**since pos does not change while we predict all 20 possible mutations for a position
* → seq_from seq_to struc_from struc_to
* stay the same,  thus we only need to copy the network outputs once
*/
       sec =copy(sec_ori);

//            float_array_2d * sec = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
//            memcpy(sec, sec_ori,1*sizeof(float_array_2d) );
//            sec->data = malloc(sizeof (float*) * sec_ori->row_size);
//            for(int i = 0; i < sec_ori->row_size; i++){
//                sec->data[i] = malloc(sizeof (float) * sec_ori->row_size);
//                memcpy(sec->data[i], sec_ori->data[i], sizeof (float) * sec_ori->col_size);
//            }
        a =copy(a_ori);
//            float_array_2d * a = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
//            memcpy(a, a_ori,1*sizeof(float_array_2d) );
//            a->data = malloc(sizeof (float*) * a_ori->row_size);
//            for(int i = 0; i < a_ori->row_size; i++){
//                a->data[i] = malloc(sizeof (float) * a_ori->row_size);
//                memcpy(a->data[i], a_ori->data[i], sizeof (float) * a_ori->col_size);
//            }

        u =copy(u_ori);
//            float_array_2d * u = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
//            memcpy(u, u_ori,1*sizeof(float_array_2d) );
//            u->data = malloc(sizeof (float*) * u_ori->row_size);
//            for(int i = 0; i < u_ori->row_size; i++){
//                u->data[i] = malloc(sizeof (float) * u_ori->row_size);
//                memcpy(u->data[i], u_ori->data[i], sizeof (float) * u_ori->col_size);
//            }

        b =copy(b_ori);
//            float_array_2d * b = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
//            memcpy(b, b_ori,1*sizeof(float_array_2d) );
//            b->data = malloc(sizeof (float*) * b_ori->row_size);
//            for(int i = 0; i < b_ori->row_size; i++){
//                b->data[i] = malloc(sizeof (float) * b_ori->row_size);
//                memcpy(b->data[i], b_ori->data[i], sizeof (float) * b_ori->col_size);
//            }


       for(int mut = 0; mut < 20; mut++){
           if(aa_ori == revers_aa[mut] && snap_20<1 )					// «mr  // allowing the switch to 20 mutations including native → native
               continue;
           char pathToSaveResult[256];
           sequ[pos] = revers_aa[mut];

           input_size=seq_to- seq_from+1;
           a_tmp=create_and_run_from_to(FA_MODEL,NULL,seq_from,seq_to);
           u_tmp=create_and_run_from_to(FU_MODEL,NULL,seq_from,seq_to);
           b_tmp=create_and_run_from_to(FB_MODEL,NULL,seq_from,seq_to);
//              a_tmp = run_model(models[FA_MODEL], create_inputs(features[FA_MODEL], seq_from, seq_to,sequ,chain_length,models[FA_MODEL]->num_input,NULL),input_size);
//              u_tmp = run_model(models[FU_MODEL], create_inputs(features[FU_MODEL], seq_from, seq_to,sequ,chain_length,models[FU_MODEL]->num_input,NULL),input_size);
//              b_tmp = run_model(models[FB_MODEL], create_inputs(features[FB_MODEL], seq_from, seq_to,sequ,chain_length,models[FB_MODEL]->num_input,NULL),input_size);

           size_t iter = 0;
           for(size_t j=seq_from;j<= seq_to;j++,iter++) {
               free(a->data[j]);
               free(u->data[j]);
               free(b->data[j]);
               a->data[j] = a_tmp->data[iter];
               u->data[j] = u_tmp->data[iter];
               b->data[j] = b_tmp->data[iter];
           }

           input_size=struc_to- struc_from+1;
           uu=create_and_run_from_to(FUU_MODEL,u,struc_from,struc_to);
           ub=create_and_run_from_to(FUB_MODEL,u,struc_from,struc_to);
           bu=create_and_run_from_to(FBU_MODEL,b,struc_from,struc_to);
           bb=create_and_run_from_to(FBB_MODEL,b,struc_from,struc_to);
//            uu = run_model(models[FUU_MODEL], create_inputs(features[FUU_MODEL], struc_from, struc_to,sequ,chain_length,models[FUU_MODEL]->num_input, u),input_size);
//             ub = run_model(models[FUB_MODEL], create_inputs(features[FUB_MODEL], struc_from, struc_to,sequ,chain_length,models[FUB_MODEL]->num_input, u),input_size);
//             bu = run_model(models[FBU_MODEL], create_inputs(features[FBU_MODEL], struc_from, struc_to,sequ,chain_length,models[FBU_MODEL]->num_input, b),input_size);
//            bb = run_model(models[FBB_MODEL], create_inputs(features[FBB_MODEL], struc_from, struc_to,sequ,chain_length,models[FBB_MODEL]->num_input, b),input_size);
            sec_tmp=jury(uu, ub, bu, bb);

           iter = 0;
           for(size_t j=struc_from;j<= struc_to;j++,iter++) {
               free(sec->data[j]);
               sec->data[j] = sec_tmp->data[iter];
           }
           //

    if(quiet<1)         printf( "%s%c%u%c         \r", "Done: ",aa_ori,pos+1,revers_aa[mut]);
	 //   printf( "%s%c%u%c\n", "Done: ",aa_ori,pos+1,revers_aa[mut]);
	   
	   
           snprintf(pathToSaveResult, sizeof pathToSaveResult, "%s_%c%u%c", pathToOut,aa_ori,pos+1,revers_aa[mut]);
	   
	   if(short_output>0){ write_short_output(pathToSaveResult, sec, a, sequ, chain_length,struc_from,struc_to );} 			//← mr: added short oouput option only the changed part+- winsize will be put to file, faster writing here and reading in snap"!
	    else           write_output(pathToSaveResult, sec, a, sequ, chain_length);
           //free run
	    
	    write_extra_out(EXTRA,  sec, a, sequ, pos );
	    
           free_float_array_2d_ptr_only(a_tmp);
           free_float_array_2d_ptr_only(u_tmp);
           free_float_array_2d_ptr_only(b_tmp);

           free_float_array_2d(uu);
           free_float_array_2d(ub);
           free_float_array_2d(bu);
           free_float_array_2d(bb);

           free_float_array_2d_ptr_only(sec_tmp);
       }
       free_float_array_2d(a);
       free_float_array_2d(u);
       free_float_array_2d(b);
       free_float_array_2d(sec);

       sequ[pos] = aa_ori;		/// restore original sequence
   }
   fclose(EXTRA);
}

 void help(){
    printf("NAME:\n");
    printf("\tmadprof\n");
    printf("\n");
    printf("DESCRIPTION:\n");
    printf("    Secondary structure and solvent accessibility prediction using neural networsk\n");
    printf("    \n");
    printf("USAGE:\n");
    printf("    \tmadprof --fasta [query.fasta] --out [path to dir] --model [path to model]\n");
    printf("    \n");
    printf("OPTIONS:\n");
    printf("    \t--fasta\n");
    printf("    \tInput (single) FASTA file\n");
    printf("    \t--out\n");
    printf("    \tEither an output file or a directory. If not provided or a directory, the suffix of the input filename is replaced to create an output filename\n");
    printf("    \t--model\n");
    printf("    \tDirectory where the model and feature files are stored\n");

    printf("EXAMPLES:\n");
    printf("\n");
    printf("    \tmadprof -fasta query.fasta \n");
    printf("    \tmadprof -fasta query.fasta -out out_dir/\n");
    printf("    \tmadprof --input query.blastPsiMat -out out_dir/\n");
    printf("    \tmadprof  -i /home/r/roosm/Mprof/893f418a5f0e2760e6041d81d058e578.blastPsiMat  --out /tmp/mprof/893f418a5f0e2760e6041d81d058e578. --model /var/tmp/mamut/Mprof_models/ \n");
    
    printf("    \tmadprof -d | --debug : give additionl dbug output \n");
    printf("    \tmadprof -s20 | --snap_20 : use 20 mutations even thow A->A probably wont make any diffrence \n"); 
    printf("    \tmadprof -s | --short : ommit unchanged portions in output faster filewrite/read/parse  \n");
       printf("    \tmadprof -q | --quiet :  none output for clusterworkers\n");
    
    exit( 0);
}

void set_snap(){snap_20=20;}

/** a lot more easy to read then the run_model( ... create inputs(....)...) thingi, just wraps*/
float_array_2d * create_and_run(int Model_num ,float_array_2d * pre_output){
	return run_model(models[Model_num], create_inputs(features[Model_num], 0, chain_length - 1,sequ,chain_length,models[Model_num]->num_input,pre_output),input_size);
}
/** a lot more easy to read then the run_model( ... create inputs(....)...) thingi, just wraps*/
float_array_2d * create_and_run_from_to(int Model_num ,float_array_2d * pre_output,size_t from,size_t to){
	return  run_model(models[Model_num], create_inputs(features[Model_num], from, to,sequ,chain_length,models[Model_num]->num_input, pre_output),input_size);
}

/** args parsing*/
	int parse_args (int argc, const char * argv[]){
	    int paramCount=0;
	    
	  for (int i=1; i<argc; i++){
	      /// switches;
	      if(argv[i][0]!='-' ) continue;
	      if(!strcmp(argv[i],"-h")|| !strcmp(argv[i],"--help")) {help();  return -1;}			// »mr // added: --help
	      if(!strcmp(argv[i],"--snap_20")|| !strcmp(argv[i],"-s20")) {set_snap();}			// »mr // added: --help
	      if(!strcmp(argv[i],"--dbug")|| !strcmp(argv[i],"-d")) {dbug++;}			// »mr // added: --help
	      if(!strcmp(argv[i],"--short")|| !strcmp(argv[i],"-s")) {short_output++;}
	      if(!strcmp(argv[i],"--quiet")|| !strcmp(argv[i],"-q")) {quiet++;}
	      /// for all other args the i may NOT reach argc -1 
	    }
	    
	    for (int i=1; i<argc-1; i++){
	      /// switches;

	      /// for all other args the i may NOT reach argc -1 
	      if(i+1 <argc){
	        if (!strcmp(argv[i],"--model")) {
	            pathToModel=argv[++i];
	        }
	        if (!strcmp(argv[i],"--out")||!strcmp(argv[i],"-o")){
	            pathToOut=argv[++i];
				// FIXME: Error handling for unsuccessful case?
		    if(!f_exists(dirname( strdup(pathToOut)))){
		      printf("Error: output file [%s] not found ... maybe try \" -o /tmp/sometmpfile \"\n",pathToOut);
		      exit( -14);
		    }
		   if(dbug>0)  printf( "-o : %s \n",pathToOut );
	        } 
	        if (!strcmp(argv[i],"--input")||!strcmp(argv[i],"-i")){
		  in_file=argv[++i];
//#in_file=get_filename(argv[i+1]);
		  if(!f_exists(in_file)){
		     printf("Error: input file [%s] not found ... maybe try \" -i /home/r/roosm/Mprof/893f418a5f0e2760e6041d81d058e578.blastPsiMat \"\n",argv[i+1]);
		     exit( -13);
		  }
		  
	            if(!strcmp(get_filename_ext(in_file),"blastPsiMat")){
	                sequ=parse_blast_pis_mat(in_file);
	                hasPsiMat=1;
	            } else {
	                /* The MIT License start */
	                sequ=parse_fasta(in_file);
	                /* The MIT License end */
	            }
	            out_file= strdup(in_file);
	            out_file=get_filename(out_file);
	            paramCount++;
	            chain_length = strlen(sequ);
		   if(dbug>0)  printf( "-i : %s \n", in_file);
	        }
	        
	//        else{	            i++;	        }	// bad habit might not add up !!
	      }
	    } // end of for-loop for command line input
	    
	    if(argc == 1 )
	    { help(); exit( 0);}
	    if(paramCount!=1){
	        printf("Error: Wrong parameter");
	        exit( 1);
	    }else{
	        if(pathToOut==NULL){
	            snprintf(outPahtBuf, sizeof outPahtBuf, "%s%s%s.reprof", pwd,"//",out_file);
	        }else{
	            snprintf(outPahtBuf, sizeof outPahtBuf, "%s", pathToOut);
	        }

	    }
		return 0;// -1;			/// we are not jet done but can start the work
	}

	/** copy things, worth a function to reduce error probability and make code readable*/
	  float_array_2d * copy(float_array_2d * _ori){
		float_array_2d * c = (float_array_2d * ) malloc(1*sizeof(float_array_2d));
		memcpy(c, _ori,1*sizeof(float_array_2d) );
		c->data = malloc(sizeof (float*) * _ori->row_size);
		for(int i = 0; i < _ori->row_size; i++){
			c->data[i] = malloc(sizeof (float) * _ori->row_size);
			memcpy(c->data[i], _ori->data[i], sizeof (float) * _ori->col_size);
		}
		return c;
	  }
	  /**function definitions*/
	  const char *get_filename_ext(const char *filepath) {
	      const char *dot = strrchr(filepath, '.');
	      if(!dot || dot == filepath) return "";
	      return dot + 1;
	  }

	   char *get_filename( char *filepath) {
	       char *filename = basename( strdup(filepath));					//← strdup added else we get some strange errs
	       char *file = strtok(filename, ".");
	      if(!file || file == filename) return filename;

	      return file;
	  }

	  features_list * parse_feature_file(const char *configuration_file)
	  {
	      feature_item * currElm, * head_list;
	      head_list = NULL;
	  	FILE *conf = fopen(configuration_file, "r");
	  	if(!conf)
	  	{
	  		return NULL;
	  	}
	      char line[200];

	      while(fgets(line,sizeof(line),conf) != NULL){
	          char head [80];
	          struct feature * feat;
	          feat = (struct feature *)malloc(sizeof(struct feature));
	          sscanf(line,"%s ",&head);
	          if(strcmp(head,"option") == 0)
	          {
	              free(feat);
	              continue;

	          }
	          if(strcmp(head,"output") == 0)
	          {
	              char windowStr [10];
	              sscanf(line, "%s %s %s %s", &head,feat->source,feat->feature,windowStr);
	              feat->window = atoi(windowStr);
	              currElm = (feature_item *)malloc(sizeof(feature_item));
	              currElm->val = feat;
	              currElm->next  = head_list;
	              head_list = currElm;
	              continue;
	          }
	          if(strcmp(head,"input") == 0)
	          {
	              char windowStr [10];
	              sscanf(line, "%s %s %s %s", &head,feat->source,feat->feature,windowStr);
	              feat->window = atoi(windowStr);
	              if (feat->window > max_window) {
	                  max_window = feat->window;
	              }
	              currElm = (feature_item *)malloc(sizeof(feature_item));
	              currElm->val = feat;
	              currElm->next  = head_list;
	              head_list = currElm;
	              continue;
	          }
	      }
	  	fclose(conf);
	      feature_item *feature=head_list;
	      currElm = NULL;
	      feature_item *new_head_list=NULL;
	      do{
	          currElm = (feature_item *)malloc(sizeof(feature_item));
	          currElm->val = (struct feature *)malloc(sizeof(struct feature));
	          memcpy(currElm->val,feature->val,sizeof(struct feature));
	          currElm->next  = new_head_list;
	          new_head_list = currElm;
	      } while((feature=feature->next)!=NULL);
	      free_feature_list(head_list);
	  	return new_head_list;
	  }
