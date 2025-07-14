import wget

def wget_download(url, dir_file, desc=None):
    """
    Downloads a file from a given URL using the wget library.

    Usage:
        >>> wget_out = wget_download('https://example.com/file.txt', '/path/to/save/file.txt', 'Downloading file.txt')

    Inputs:
        url -> [str] URL of the file to be downloaded.
        dir_file -> [str] Path where the downloaded file will be stored.
        desc -> [str, optional, default=None] Description of the download. If provided, it will be printed before the download starts.
    Outputs:
        wget_out -> [str] Path of the downloaded file.
    """

    # Print the provided description
    if desc: print(desc)

    # Use wget to download the file and store it at the specified path
    wget_out = wget.download(url, dir_file)

    return wget_out


