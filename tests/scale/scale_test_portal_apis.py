import random
import time

from locust import HttpUser, between, task

random.seed(time.time())

DEV_COLLECTION_IDS = ["673637cf-dcb7-45e1-bb88-72a27c50c8ca",
                      "a55cac37-629f-4d28-9b3f-4ece7a7ac174",
                      "e7c47289-22fe-440e-ba24-a4434ac3f379",
                      "18e6337b-8c36-4977-9ac8-ca8aa5494481"]


class WebsiteUser(HttpUser):
    wait_time = between(1, 2)

    @task
    def get_collections(self):
        self.client.get("dp/v1/collections")

    @task
    def get_collection_info(self):
        self.client.get(f"dp/v1/collections/{random.choice(DEV_COLLECTION_IDS)}")
