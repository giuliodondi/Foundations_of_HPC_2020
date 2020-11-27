#!/bin/bash

LC_ALL='en_US.UTF-8'

#prepare the usage printout to display with the -h or --help  option or error messages
print_userguide()
{
	echo "Generates random integer numbers."
	echo "Mandatory arguments:"
	echo "-n, --num : Specifies the amount of numbers to generate. Must be a positive integer."
	echo "Optional Arguments:"
	echo "-s, --seed : Specifies the seed to use for generation. Must be a positive integer. Defaults to the process ID"
	echo "-m, --max  : Specifies the maximum value for the random numbers. Must be a positive integer." 
	echo "-f, --float: Generates random decimal numbers"
}


throw_missingarg_err()
{
	case $1 in
	1)
	echo "ERROR: missing a required argument."
	;;
	2)
	echo "ERROR: unspecified or invalid option."
	;;
	*);;
	esac

	echo "USe the -h or --help options for information."
	exit $1

}

#TODO: handle the case where the variable is not a number
manage_arg()
{

	if [ -z ${1} ]
	then 
	throw_missingarg_err 1
	elif [ ${1} -lt 1 ]
	then
	throw_missingarg_err 2
	fi

}

rand_float()
{
	${1}
	leftdigits=${new_num}
	${2}
	rightdigits=${new_num}
	new_num=${leftdigits}'.'${rightdigits}
}

rand_unbound()
{
	new_num=${RANDOM}
}

rand_bound()
{
	new_num=$((RANDOM%VAR_M))
}




#check for missing required arguments
if [[ $# -eq 0 ]]
then
        throw_missingarg_err 1
else
while [[ $# -gt 0 ]]
do
	case $1 in 
		-h|--help)
			print_userguide
			exit 0
		;;
		-n|--num)
			shift
			VAR_N=${1}
			manage_arg ${VAR_N}
		;;
		-s|--seed)
			shift
			VAR_S=${1}
			manage_arg ${VAR_S}
		;;
		-m|--max)
			shift
			VAR_M=${1}
			manage_arg ${VAR_M}
		;;
		-f|--float)
			VAR_F='true'
		;;
		*)
			throw_missingarg_err 2
	esac
	shift	#so that $1 parses all the arguments one by one
done
fi

#check for missing required arguments
if [ -z ${VAR_N} ]
then
        throw_missingarg_err 1
fi

#set the seed
if [ -z ${VAR_S} ]
then
	RANDOM=$$
else
	RANDOM=${VAR_S}
fi


#define the generator function to call
#based on whether the bound is specified and whether the option is int or float 

if [ -z ${VAR_F} ]
then
	#integer number generation
	if [ -z ${VAR_M} ]
	then
		GENFUN="rand_unbound"

	else
		GENFUN="rand_bound"
	fi
	formatspec="%d: %d\n"
else
	#float number generation
	if [ -z ${VAR_M} ]
	then
		f1="rand_unbound"

	else
		let 'VAR_M=--VAR_M'
		f1="rand_bound"
	fi	
	f2="rand_unbound"
	GENFUN="rand_float ${f1} ${f2}"
	formatspec="%d: %f\n"
fi



#generate the random numbers
printf "Generating %s random integer(s)\n" "${VAR_N}"

new_num=0

for (( i=1; i<=${VAR_N}; i++ ))
do

	${GENFUN}

	printf "${formatspec}" ${i} ${new_num}
done

