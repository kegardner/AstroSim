import sqlite3 as sql

def connect():
    try:
        # connecting to the sq;ite db
        with sql.connect("..\data") as conn:
            print('Connected to the SQLite DB.')
            return conn
    except (sql.DatabaseError, Exception) as error:
        print(error)

if __name__ == '__main__':
    conn = connect()
    cur = conn.cursor()

    cur.execute("select * from fb_redshift;")
    result = cur.fetchall()
    print(result)