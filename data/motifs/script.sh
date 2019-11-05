#!/bin/bash
# param -h for HOCOMOCO, -j for JASPAR, -all for all (also mtf for TRANSFAC and UniPROBE are created)

# check parameters

if [ $# -eq 0 ]
    then
         read -p "Type -h to use HOCOMOCO, -j to use JASPAR and -all to use all " source
elif [ $# -eq 1 ]
    then
        if [ $1 == "-h" ] || [ $1 == "-j" ] || [ $1 == "-all" ]
            then
                source=$1
        else
            read -rp "Type -h to use HOCOMOCO, -j to use JASPAR and -all to use all " source
        fi
elif [ $# -eq 2 ]
    then
        if [ $1 == "-h" ] && [ $2 == "-j" ]
            then
                source="-all"
        elif [ $1 == "-j" ] && [ $2 == "-h" ]
            then
                source="-all"
        else
           read -rp "Type -h to use HOCOMOCO, -j to use JASPAR and -all to use all " source
        fi
else
    read -rp "Type -h to use HOCOMOCO, -j to use JASPAR and -all to use all " source
fi

# if $source is set correctly, start creating all desired files(Pwm, Mtf)

if [ "${source}" == "-h" ]
    then
        echo
        echo "HOCOMOCO"
        echo
        echo "Create Pwm files..."
        wget -c http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_pcms_HUMAN_mono.txt
        wget -c http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_pcms_MOUSE_mono.txt
        python createPwm.py -i HOCOMOCOv11_full_pcms_HUMAN_mono.txt -f hocomoco-pcm -o hocomoco
        python createPwm.py -i HOCOMOCOv11_full_pcms_MOUSE_mono.txt -f hocomoco-pcm -o hocomoco
        echo
        echo "Create Mtf file..."
        ./make_annotation.sh -h
        echo

        # clean up
        rm HOCOMOCOv11_full_pcms_HUMAN_mono.txt
        rm HOCOMOCOv11_full_pcms_MOUSE_mono.txt



elif [ "${source}" == "-j" ]
    then
        echo
        echo "JASPAR"
        echo
        echo "Create Pwm files.."
        wget -c --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6" http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt
        wget -c --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6" http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_plants_non-redundant_pfms_jaspar.txt
        python createPwm.py -i JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt -f jaspar-2016 -o jaspar_vertebrates
        python createPwm.py -i JASPAR2020_CORE_plants_non-redundant_pfms_jaspar.txt -f jaspar-2016 -o jaspar_plants
        echo
        echo "Create Mtf file..."
        echo
        ./make_annotation.sh -j
        echo

        # clean up
        rm JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt
        rm JASPAR2020_CORE_plants_non-redundant_pfms_jaspar.txt

elif [ "${source}" == "-all" ]
    then
        echo
        echo "HOCOMOCO"
        echo
        echo "Create Pwm files..."
        wget -c http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_pcms_HUMAN_mono.txt
        wget -c http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_pcms_MOUSE_mono.txt
        python createPwm.py -i HOCOMOCOv11_full_pcms_HUMAN_mono.txt -f hocomoco-pcm -o hocomoco
        python createPwm.py -i HOCOMOCOv11_full_pcms_MOUSE_mono.txt -f hocomoco-pcm -o hocomoco
        echo
        echo "Create Mtf file..."
        echo
        ./make_annotation.sh -h
        echo

        # clean up
        rm HOCOMOCOv11_full_pcms_HUMAN_mono.txt
        rm HOCOMOCOv11_full_pcms_MOUSE_mono.txt

        echo
        echo "JASPAR"
        echo
        echo "Create Pwm files.."
        wget -c --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6" http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt
        wget -c --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6" http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_plants_non-redundant_pfms_jaspar.txt
        python createPwm.py -i JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt -f jaspar-2016 -o jaspar_vertebrates
        python createPwm.py -i JASPAR2020_CORE_plants_non-redundant_pfms_jaspar.txt -f jaspar-2016 -o jaspar_plants
        echo
        echo "Create Mtf file..."
        ./make_annotation.sh -j
        echo

        # clean up
        rm JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt
        rm JASPAR2020_CORE_plants_non-redundant_pfms_jaspar.txt

        echo
        echo "TRANSFAC"
        echo
        echo "Create Mtf file..."
        python createMtf.py --t
        echo

        echo
        echo "UniProbe"
        echo
        echo "Create Mtf file..."
        python createMtf.py --up --us
        echo

else
    echo "ERROR"
fi
