from datetime import datetime, timedelta
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from airflow import DAG
from airflow.operators.python import PythonOperator
from airflow.providers.postgres.hooks.postgres import PostgresHook
from airflow.providers.amazon.aws.hooks.s3 import S3Hook
from airflow.models import Variable

POSTGRES_CONN_ID = Variable.get("POSTGRES_CONN_ID")
AWS_CONN_ID = Variable.get("AWS_CONN_ID")
BUCKET_NAME = Variable.get("bucket_name")

def extract_data(ti, **kwargs):
    query = """
    SELECT smiles, other_columns
    FROM your_table
    WHERE date_column = '{{ ds }}';
    """
    pg_hook = PostgresHook(postgres_conn_id=POSTGRES_CONN_ID)
    df = pg_hook.get_pandas_df(sql=query)
    ti.xcom_push(key='extracted_data', value=df)

def transform_data(ti, **kwargs):
    df = ti.xcom_pull(key='extracted_data', task_ids='extract_data')
    df['mol'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(x))
    df['Molecular_Weight'] = df['mol'].apply(lambda x: Descriptors.MolWt(x))
    df['LogP'] = df['mol'].apply(lambda x: Descriptors.MolLogP(x))
    df['TPSA'] = df['mol'].apply(lambda x: Descriptors.TPSA(x))
    df['H_Donors'] = df['mol'].apply(lambda x: Descriptors.NumHDonors(x))
    df['H_Acceptors'] = df['mol'].apply(lambda x: Descriptors.NumHAcceptors(x))
    df['Lipinski_pass'] = df.apply(lambda row: (row['H_Donors'] <= 5) and (row['H_Acceptors'] <= 10) and (row['LogP'] <= 5) and (row['Molecular_Weight'] <= 500), axis=1)
    df.drop(columns=['mol'], inplace=True)  # Remove the RDKit Mol objects
    ti.xcom_push(key='transformed_data', value=df)

def save_and_load_to_s3(ti, **kwargs):
    df = ti.xcom_pull(key='transformed_data', task_ids='transform_data')
    filepath = '/tmp/processed_data.xlsx'
    df.to_excel(filepath, index=False, engine='openpyxl')
    s3_hook = S3Hook(aws_conn_id=AWS_CONN_ID)
    s3_hook.load_file(filename=filepath, key='data/processed_data.xlsx', bucket_name=BUCKET_NAME, replace=True)

default_args = {
    'owner': 'airflow',
    'depends_on_past': False,
    'email_on_failure': False,
    'email_on_retry': False,
    'retries': 1,
    'retry_delay': timedelta(minutes=5),
}

with DAG(
    'molecule_data_processing',
    default_args=default_args,
    description='A DAG for molecule data processing',
    schedule_interval=timedelta(days=1),
    start_date=datetime(2023, 1, 1),
    catchup=False,
) as dag:

    extract = PythonOperator(
        task_id='extract_data',
        python_callable=extract_data
    )

    transform = PythonOperator(
        task_id='transform_data',
        python_callable=transform_data
    )

    save_load = PythonOperator(
        task_id='save_and_load_to_s3',
        python_callable=save_and_load_to_s3
    )

    extract >> transform >> save_load
