from backend.common.server.config import create_api_app

app = create_api_app(
    api_paths_and_spec_files=[("/de", "de/api/de-api.yml")],
)

if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=True)
