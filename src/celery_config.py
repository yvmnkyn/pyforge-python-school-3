from celery import Celery

# Create a Celery instance
celery_app = Celery(
    'celery_app',
    broker='redis://localhost:6379/0',
    backend='redis://localhost:6379/0'
)

celery_app.conf.task_routes = {
    'molecules.tasks.substructure_search_task': {'queue': 'substructure_search_queue'}
}

celery_app.conf.update(
    task_serializer='json',
    accept_content=['json'],
    result_serializer='json',
    timezone='UTC',
    enable_utc=True,
)
