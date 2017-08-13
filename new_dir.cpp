#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstring>


#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

int main ()
{
    char n[20];

    unsigned char isFolder = 0x4;
     unsigned char isFile = 0x8;
    DIR *dir;
    struct dirent *ent;
    dir = opendir ("/home");
        if (dir != NULL) {

          /* print all the files and directories within directory */
          while ((ent = readdir (dir)) != NULL) {
            //printf ("%s\n", ent->d_name);

            //folder sign
            if(ent->d_type == isFolder)
            {
                cout <<ent->d_name <<"\n";


            }


          }
          closedir (dir);


        } else {
          /* could not open directory */
          perror ("");
          return 0;
        }

        cout << "=========" << endl;

}