/* Mock response from /dataset-metadata endpoint. */
export const DATASET_METADATA_RESPONSE = {
  metadata: {
    collection_contact_email: "firstlast@test.com",
    collection_contact_name: "First Last",
    collection_datasets: [
      {
        cell_count: 16245,
        collection_id: "b3d9199e-a6a6-4faf-b456-d109a3003a40",
        dataset_deployments: [
          {
            url: "https://test.com/e/test0.cxg/",
          },
        ],
        id: "498a9066-58dc-4e3a-8bfa-0d4331ea918e",
        name: "Aliquam eu porttitor enim, sit amet blandit nulla",
      },
      {
        cell_count: 24245,
        collection_id: "648bf25c-6d3f-4923-ab7b-fbcfc5d65776",
        dataset_deployments: [
          {
            url: "https://test.com/e/test1.cxg/",
          },
        ],
        id: "c9897529-dc61-42ca-af89-7cf10e7f9267",
        name: "Sed eu nisi condimentum",
      },
    ],
    collection_description:
      "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aenean quis lectus ut felis ornare lacinia vel nec erat.",
    collection_links: [
      {
        link_name: "GEO",
        link_type: "RAW_DATA",
        link_url: "https://test.com",
      },
      {
        link_name: "Synapse",
        link_type: "OTHER",
        link_url: "https://test.com",
      },
    ],
    collection_name: "Lorem ipsum dolor sit amet",
    collection_url: "https://test.com",
    dataset_id: "498a9066-58dc-4e3a-8bfa-0d4331ea918e",
    dataset_name: "Nullam ultrices urna nec congue aliquam",
  },
};
