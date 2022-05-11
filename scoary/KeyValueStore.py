import os
import logging
import sqlite3


class KeyValueStore:
    table_name: str

    def __init__(self, table_name, db_path: str = None):
        self.table_name = table_name

        if db_path is None:
            if 'KEY_VALUE_STORE_DB' in os.environ:
                db_path = os.environ['KEY_VALUE_STORE_DB']
            else:
                db_path = os.path.expanduser('~/.cache/keyvaluestore.db')

        self._db_path = db_path
        self.con, self.cur = self.get_cur()
        self.create_db()

    def __str__(self):
        return f'KeyValueStore {self.table_name} ({self._db_path})'

    def get_cur(self):
        try:
            con = sqlite3.connect(self._db_path)
            cur = con.cursor()
        except Exception as e:
            logging.warning(f'Failed to connect to db: {self._db_path}')
            raise e
        return con, cur

    def __del__(self):
        try:
            self.cur.close()
            self.con.close()
        except Exception:
            pass

    def create_db(self):
        raise NotImplementedError(f'Users of the abstract class {self.__class__} must implement this function!')

    @staticmethod
    def list_to_string(l) -> str:
        return ', '.join(f"'{e}'" for e in l)

    @staticmethod
    def list_to_string_bracket(l):
        return ', '.join(f"('{e}')" for e in l)

    def _create_db(self, columns: {str: str}, pk_col: str):
        columns = ', '.join(f'{col_name} {col_type}' for col_name, col_type in columns.items())
        sql = f'''
            CREATE TABLE IF NOT EXISTS {self.table_name} (
                {columns},
                PRIMARY KEY ({pk_col})
            );
        '''
        try:
            self.cur.execute(sql)
        except sqlite3.OperationalError as e:
            logging.warning(f'Failed to run this SQL command on db {self._db_path}:\n{sql}')
            raise e

    def drop_db(self):
        self.cur.execute(f'''DROP TABLE {self.table_name}''')
