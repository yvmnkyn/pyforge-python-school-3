import boto3
import json
import os

session = boto3.Session(
    aws_access_key_id=os.environ['AWS_ACCESS_KEY_ID'],
    aws_secret_access_key=os.environ['AWS_SECRET_ACCESS_KEY'],
    region_name='eu-north-1'
)

lambda_client = session.client('lambda')

lambda_client = session.client('lambda')

# Defining the event to send to Lambda function
event = {
    "name": "Alice"
}

# Invoke the Lambda function
response = lambda_client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event)
)

response_payload = response['Payload'].read()
print(json.loads(response_payload))
