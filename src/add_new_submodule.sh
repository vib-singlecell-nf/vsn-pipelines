#!/bin/bash

function usage {
    echo "usage: $0 repository-url [-cdh] [-b branch/tag]"
    echo "  -b branch/tag    link to branch/tag"
    echo "  -c               add module and generate commit"
    echo "  -d               link to develop branch"
    echo "  -h               display help"
    exit 1
}

if [[ $# -eq 0 ]]; then
    usage
fi

if [ `basename $PWD` != 'src' ]; then   
    echo "ERROR: You must be in the vsn-pipelines src directory to run this script"
    exit 1
fi

unset BRANCH
unset COMMIT
declare -a ARGS

while [ $# -gt 0 ]
do
    unset OPTIND
    unset OPTARG
    while getopts ':b:cdh' c
    do
    case $c in
            b)  if [[ ! -z "$BRANCH" ]]; then
                    echo "ERROR: Cannot specify -b and -d"
                    exit 1
                fi
                BRANCH=$OPTARG
                ;;
            c)  COMMIT=true
                ;;
            d)  if [[ ! -z "$BRANCH" ]]; then
                    echo "ERROR: Cannot specify -b and -d"
                    exit 1
                fi
                BRANCH=develop
                ;;
            h)  usage
                ;;
            :)  echo "$0: -$OPTARG needs a value" >&2;
                exit 1
                ;;
            \?) echo "$0: unknown option -$OPTARG" >&2;
                exit 1
                ;;
        esac
    done
    shift $((OPTIND-1))
    ARGS+=($1)
    shift
done

URL=${ARGS[0]}
REPO_NAME=`echo ${URL:0:-4} | sed 's!.*/!!'`

if [[ $URL =~ http.*:\/\/github.com\/vib-singlecell-nf\/.*\.git ]]; then
    echo "WARNING: ${URL} is a http(s) github repository! Using the the following SSH URL..."
    URL=`echo ${URL} | sed 's!http.*github.com/!git@github.com:!'`
    echo "    ${URL}"
fi

if [[ ! $URL =~ git@github.com:vib-singlecell-nf\/.*\.git ]]; then
    echo "ERROR: ${URL} is not a valid SSH address for a vib-singlecell-nf github repository! "
    exit 1
fi

if [[ ! -z "$BRANCH" ]]; then
    echo "Adding requested submodule on ${BRANCH} branch..."
    git submodule add --branch $BRANCH $URL
else
    echo "Adding requested submodule..."
    git submodule add $URL
fi

echo "Updating submodules..."
git submodule update --init --recursive

if [[ $COMMIT == 'true' ]]; then
    echo '-c passed. Adding module and commiting.'
    git add ./${REPO_NAME}
    git commit -m "Add ${REPO_NAME} submodule"
fi