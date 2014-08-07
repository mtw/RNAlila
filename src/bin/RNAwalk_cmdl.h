/** @file RNAwalk_cmdl.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.5
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef RNAWALK_CMDL_H
#define RNAWALK_CMDL_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef RNAWALK_CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define RNAWALK_CMDLINE_PARSER_PACKAGE "RNAwalk"
#endif

#ifndef RNAWALK_CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define RNAWALK_CMDLINE_PARSER_PACKAGE_NAME "RNAwalk"
#endif

#ifndef RNAWALK_CMDLINE_PARSER_VERSION
/** @brief the program version */
#define RNAWALK_CMDLINE_PARSER_VERSION VERSION
#endif

/** @brief Where the command line options are stored */
struct RNAwalk_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *full_help_help; /**< @brief Print help, including hidden options, and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * walktype_arg;	/**< @brief Specify type of walk.
  A ... Adaptive walk
  G ... Gradient Walk
  R ... Random walk
 (default='G').  */
  char * walktype_orig;	/**< @brief Specify type of walk.
  A ... Adaptive walk
  G ... Gradient Walk
  R ... Random walk
 original value given at command line.  */
  const char *walktype_help; /**< @brief Specify type of walk.
  A ... Adaptive walk
  G ... Gradient Walk
  R ... Random walk
 help description.  */
  int walklength_arg;	/**< @brief Specify walk length. If the walk ends before the predefined numberof steps is reached (e.g. in case a local minimum has been reached inan adaptive or gradient walk) is reached, the walk will terminate.
 (default='100000').  */
  char * walklength_orig;	/**< @brief Specify walk length. If the walk ends before the predefined numberof steps is reached (e.g. in case a local minimum has been reached inan adaptive or gradient walk) is reached, the walk will terminate.
 original value given at command line.  */
  const char *walklength_help; /**< @brief Specify walk length. If the walk ends before the predefined numberof steps is reached (e.g. in case a local minimum has been reached inan adaptive or gradient walk) is reached, the walk will terminate.
 help description.  */
  double temp_arg;	/**< @brief Rescale energy parameters to a temperature of temp C.
 (default='37').  */
  char * temp_orig;	/**< @brief Rescale energy parameters to a temperature of temp C.
 original value given at command line.  */
  const char *temp_help; /**< @brief Rescale energy parameters to a temperature of temp C.
 help description.  */
  int dangles_arg;	/**< @brief How to treat \"dangling end\" energies for bases adjacent to helicesin free ends and multi-loops
 (default='2').  */
  char * dangles_orig;	/**< @brief How to treat \"dangling end\" energies for bases adjacent to helicesin free ends and multi-loops
 original value given at command line.  */
  const char *dangles_help; /**< @brief How to treat \"dangling end\" energies for bases adjacent to helicesin free ends and multi-loops
 help description.  */
  int noLP_flag;	/**< @brief Produce structures without lonely pairs (helices of length 1).
 (default=off).  */
  const char *noLP_help; /**< @brief Produce structures without lonely pairs (helices of length 1).
 help description.  */
  double betaScale_arg;	/**< @brief Set the scaling of the Boltzmann factors.
 (default='1.').  */
  char * betaScale_orig;	/**< @brief Set the scaling of the Boltzmann factors.
 original value given at command line.  */
  const char *betaScale_help; /**< @brief Set the scaling of the Boltzmann factors.
 help description.  */
  int verbose_flag;	/**< @brief Verbose output
 (default=off).  */
  const char *verbose_help; /**< @brief Verbose output
 help description.  */
  int debug_flag;	/**< @brief Debugging output
 (default=off).  */
  const char *debug_help; /**< @brief Debugging output
 help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int full_help_given ;	/**< @brief Whether full-help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int walktype_given ;	/**< @brief Whether walktype was given.  */
  unsigned int walklength_given ;	/**< @brief Whether walklength was given.  */
  unsigned int temp_given ;	/**< @brief Whether temp was given.  */
  unsigned int dangles_given ;	/**< @brief Whether dangles was given.  */
  unsigned int noLP_given ;	/**< @brief Whether noLP was given.  */
  unsigned int betaScale_given ;	/**< @brief Whether betaScale was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int debug_given ;	/**< @brief Whether debug was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct RNAwalk_cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure RNAwalk_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure RNAwalk_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *RNAwalk_args_info_purpose;
/** @brief the usage string of the program */
extern const char *RNAwalk_args_info_usage;
/** @brief all the lines making the help output */
extern const char *RNAwalk_args_info_help[];
/** @brief all the lines making the full help output (including hidden options) */
extern const char *RNAwalk_args_info_full_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int RNAwalk_cmdline_parser (int argc, char **argv,
  struct RNAwalk_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use RNAwalk_cmdline_parser_ext() instead
 */
int RNAwalk_cmdline_parser2 (int argc, char **argv,
  struct RNAwalk_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int RNAwalk_cmdline_parser_ext (int argc, char **argv,
  struct RNAwalk_args_info *args_info,
  struct RNAwalk_cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int RNAwalk_cmdline_parser_dump(FILE *outfile,
  struct RNAwalk_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int RNAwalk_cmdline_parser_file_save(const char *filename,
  struct RNAwalk_args_info *args_info);

/**
 * Print the help
 */
void RNAwalk_cmdline_parser_print_help(void);
/**
 * Print the full help (including hidden options)
 */
void RNAwalk_cmdline_parser_print_full_help(void);
/**
 * Print the version
 */
void RNAwalk_cmdline_parser_print_version(void);

/**
 * Initializes all the fields a RNAwalk_cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void RNAwalk_cmdline_parser_params_init(struct RNAwalk_cmdline_parser_params *params);

/**
 * Allocates dynamically a RNAwalk_cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized RNAwalk_cmdline_parser_params structure
 */
struct RNAwalk_cmdline_parser_params *RNAwalk_cmdline_parser_params_create(void);

/**
 * Initializes the passed RNAwalk_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void RNAwalk_cmdline_parser_init (struct RNAwalk_args_info *args_info);
/**
 * Deallocates the string fields of the RNAwalk_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void RNAwalk_cmdline_parser_free (struct RNAwalk_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int RNAwalk_cmdline_parser_required (struct RNAwalk_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* RNAWALK_CMDL_H */
