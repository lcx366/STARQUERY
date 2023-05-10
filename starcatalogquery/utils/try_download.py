import wget

def wget_download(url,dir_file,desc=None):
    """
    Download a file using the wget command line.

    Inputs:
        url -> [str] URL of the file to download
        dir_file -> [str] Path of the file downloaded
        desc -> [str,optional,default=None] Description of the downloading   

    Outputs:
        wget_out -> [str] Same as the dir_file 
    """
    if desc: print(desc)
    wget_out = wget.download(url,dir_file)
    print()

    return wget_out