import wc_utils.quilt


def main():
    '''Backup or download data from/to Quilt
    '''
    path = input("BSON file location:\n")
    package = 'datanator'
    manager = wc_utils.quilt.QuiltManager(path=path, package=package)
    backup = input("Backup or Download (choose 'backup' or 'download')?\n")    
    if backup.lower() == 'backup':
        message = input("Optionally, enter a commit message:\n")
        manager.upload_package(message=message or None)
    else:
        manager.download_package()


if __name__ == '__main__':
    main()
