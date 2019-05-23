import wc_utils.quilt


def main():
    '''Backup or download data from/to Quiltdata
    '''
    path = input("BSON file location: \n")
    pakcage = 'datanator_nosql'
    token = 'eyJpZCI6ICJiNGZkOWUzZS03ZWIyLTQ2NzAtYjdkNy02MTgyNjE2ZWI4ZTIiLCAiY29kZSI6ICI1MTUwZDU2Zi1iYTY2LTRmZWQtODM5Yi0zMzg5NTBiN2Q1ZTUifQ=='
    manager = wc_utils.quilt.QuiltManager(path, pakcage, token=token,verbose=True)
    print(manager.token)
    backup = input("Backup or Download (choose 'backup' or 'download')?\n")
    if backup.lower() == 'backup':
        manager.upload()
    else:
        manager.download()


if __name__ == '__main__':
    main()
