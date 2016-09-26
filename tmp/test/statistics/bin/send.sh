sendEmail   -f  catrapid.crg@gmail.com        \
              -t  $1                  \
              -s  smtp.gmail.com:587  \
              -xu catrapid.crg@gmail.com        \
              -xp gianga1976          \
              -u "catRAPID results with seed"  $2      \
              -m "Hi Buddy,\n \n Your calculations are ready at: http://"$3":8888/tmp/download/report."$2".pdf \n \n The catRAPID team. \n \n ---\n Two posts available in our group: 1 Experimental post-doctoral researcher and 1 Laboratory technician. Both positions are intended for investigating protein-RNA interactions. Please send a message to catrapid.crg@gmail.com if interested."
