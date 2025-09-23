# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is the CZ CELLxGENE Discover data portal - a platform for publishing, discovering and exploring single-cell datasets. It's a full-stack application with Python backend services and a Next.js frontend.

## Development Commands

### Python/Backend Commands

- `make fmt` - Auto-format Python code with black and run pre-commit hooks
- `make lint` - Run ruff linting on Python code
- `make unit-test` - Run all unit tests (equivalent to `make local-unit-test`)
- `make functional-test` - Run functional tests with pytest
- `make local-init` - Set up complete local development environment with Docker
- `make local-start` - Start existing local Docker environment
- `make local-stop` - Stop local Docker environment
- `make local-logs CONTAINER=<name>` - View logs for specific container
- `make local-shell CONTAINER=<name>` - Open shell in container
- `make local-dbconsole` - Connect to local PostgreSQL database

### Frontend Commands (run from `/frontend`)

- `npm run dev` - Start development server with HTTPS
- `npm run build` - Build production bundle
- `npm run start` - Start production server
- `npm run lint` - Run ESLint, Stylelint, and TypeScript checks
- `npm run type-check` - Run TypeScript type checking only
- `npm run format` - Auto-format code with Prettier and ESLint --fix
- `npm run e2e` - Run Playwright e2e tests
- `npm test` - Run Playwright tests (alias for e2e)

### Environment Setup

- Set `AWS_PROFILE=single-cell-dev` for local development
- Set `DEPLOYMENT_STAGE=test` for local Docker environment
- Copy `frontend/src/configs/local.js` to `frontend/src/configs/configs.js` for frontend development
- Run `pre-commit install` to set up git hooks

## Architecture

### Backend Structure

- **`/backend/api_server/`** - Main API server with Flask/Connexion
- **`/backend/layers/`** - Shared utilities and database models
- **`/backend/common/`** - Common utilities across services
- **`/backend/wmg/`** - Where's My Gene functionality
- **`/backend/de/`** - Differential expression analysis
- **`/backend/cellguide/`** - Cell Guide pipeline and services
- **`/backend/database/`** - Database schema and migrations
- **`/backend/portal/`** - Core portal business logic

### Frontend Structure

- **Next.js application** with TypeScript and React 18
- **Styling**: Emotion CSS-in-JS, Material-UI components, Blueprint.js, CZI SDS components
- **Key directories**:
  - `/src/pages/` - Next.js pages and routing
  - `/src/views/` - Main view components
  - `/src/components/` - Shared React components
  - `/src/common/` - Utilities and shared logic

### Development Environments

1. **Local Python-only**: Fast testing environment without Docker
2. **Local Docker**: Full environment with all services (recommended for integration work)
3. **Remote deployment (rDev)**: AWS-hosted development environment

### Database & Testing

- **PostgreSQL** for primary data storage
- **Unit tests**: Run in Docker containers with coverage reporting
- **Functional tests**: Test against deployed environments
- **E2E tests**: Playwright-based frontend testing

### Key Technologies

- **Backend**: Python 3.10, Flask, SQLAlchemy, Alembic, pytest
- **Frontend**: Next.js 14, React 18, TypeScript, Material-UI, Emotion
- **Infrastructure**: Docker, AWS, PostgreSQL
- **CI/CD**: GitHub Actions with pre-commit hooks (black, ruff, prettier)

### Code Quality

- Python: black formatting, ruff linting
- Frontend: Prettier formatting, ESLint linting, TypeScript checking
- Pre-commit hooks enforce formatting standards
- All tests must pass before deployment
