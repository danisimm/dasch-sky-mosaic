# Container image for the DASCH -> WWT link Lambda.
# AWS's official Lambda base image: Amazon Linux + Python 3.12 + the Lambda runtime.
FROM public.ecr.aws/lambda/python:3.12

COPY requirements.txt ${LAMBDA_TASK_ROOT}/
RUN pip install --no-cache-dir -r ${LAMBDA_TASK_ROOT}/requirements.txt

COPY lambda_function.py ${LAMBDA_TASK_ROOT}/

# file.function that Lambda calls for each request
CMD ["lambda_function.lambda_handler"]
