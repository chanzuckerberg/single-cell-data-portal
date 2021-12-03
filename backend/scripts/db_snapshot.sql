-- Generated via: pg_dump -Fp --no-privileges --no-owner --no-comments --schema-only
--
-- PostgreSQL database dump
--

-- Dumped from database version 10.11
-- Dumped by pg_dump version 11.9

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: collectionvisibility; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.collectionvisibility AS ENUM (
    'PRIVATE',
    'PUBLIC'
);


--
-- Name: conversionstatus; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.conversionstatus AS ENUM (
    'CONVERTED',
    'CONVERTING',
    'FAILED',
    'NA'
);


--
-- Name: datasetartifactfiletype; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.datasetartifactfiletype AS ENUM (
    'H5AD',
    'RDS',
    'CXG'
);


--
-- Name: datasetartifacttype; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.datasetartifacttype AS ENUM (
    'ORIGINAL',
    'REMIX'
);


--
-- Name: projectlink; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.projectlink AS ENUM (
    'DOI',
    'LAB_WEBSITE',
    'OTHER',
    'PROTOCOL',
    'RAW_DATA'
);


--
-- Name: uploadstatus; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.uploadstatus AS ENUM (
    'CANCELED',
    'CANCEL_PENDING',
    'FAILED',
    'NA',
    'UPLOADED',
    'UPLOADING',
    'WAITING'
);


--
-- Name: validationstatus; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.validationstatus AS ENUM (
    'INVALID',
    'NA',
    'VALID',
    'VALIDATING'
);


SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: alembic_version; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.alembic_version (
    version_num character varying(32) NOT NULL
);

--
-- Data for Name: alembic_version; Type: TABLE DATA; Schema: public; Owner: corpora_dev
--

COPY public.alembic_version (version_num) FROM stdin;
09_7794b1ea430f
\.

--
-- Name: dataset; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.dataset (
    id character varying NOT NULL,
    revision integer,
    name character varying,
    organism jsonb,
    tissue jsonb,
    assay jsonb,
    disease jsonb,
    sex jsonb,
    ethnicity jsonb,
    created_at timestamp without time zone DEFAULT now() NOT NULL,
    updated_at timestamp without time zone DEFAULT now() NOT NULL,
    development_stage jsonb,
    collection_id character varying NOT NULL,
    collection_visibility public.collectionvisibility NOT NULL,
    cell_count integer,
    is_valid boolean
);


--
-- Name: dataset_artifact; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.dataset_artifact (
    id character varying NOT NULL,
    dataset_id character varying NOT NULL,
    filename character varying,
    filetype public.datasetartifactfiletype,
    type public.datasetartifacttype,
    user_submitted boolean,
    s3_uri character varying,
    created_at timestamp without time zone DEFAULT now() NOT NULL,
    updated_at timestamp without time zone DEFAULT now() NOT NULL
);


--
-- Name: dataset_processing_status; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.dataset_processing_status (
    id character varying NOT NULL,
    dataset_id character varying NOT NULL,
    upload_status public.uploadstatus,
    upload_progress double precision,
    upload_message character varying,
    validation_status public.validationstatus,
    validation_message character varying,
    rds_status public.conversionstatus,
    cxg_status public.conversionstatus,
    h5ad_status public.conversionstatus,
    created_at timestamp without time zone DEFAULT now() NOT NULL,
    updated_at timestamp without time zone DEFAULT now() NOT NULL
);



--
-- Name: project; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.project (
    id character varying NOT NULL,
    owner character varying NOT NULL,
    name character varying,
    description character varying,
    created_at timestamp without time zone DEFAULT now() NOT NULL,
    updated_at timestamp without time zone DEFAULT now() NOT NULL,
    visibility public.collectionvisibility NOT NULL,
    contact_email character varying,
    contact_name character varying,
    data_submission_policy_version character varying NOT NULL,
    obfuscated_uuid character varying
);


--
-- Name: project_link; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.project_link (
    id character varying NOT NULL,
    collection_id character varying NOT NULL,
    link_url character varying,
    link_type public.projectlink,
    created_at timestamp without time zone DEFAULT now() NOT NULL,
    updated_at timestamp without time zone DEFAULT now() NOT NULL,
    link_name character varying,
    collection_visibility public.collectionvisibility NOT NULL
);


--
-- Name: alembic_version alembic_version_pkc; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.alembic_version
    ADD CONSTRAINT alembic_version_pkc PRIMARY KEY (version_num);


--
-- Name: dataset_artifact dataset_artifact_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.dataset_artifact
    ADD CONSTRAINT dataset_artifact_pkey PRIMARY KEY (id);


--
-- Name: dataset dataset_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.dataset
    ADD CONSTRAINT dataset_pkey PRIMARY KEY (id);


--
-- Name: dataset_processing_status dataset_processing_status_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.dataset_processing_status
    ADD CONSTRAINT dataset_processing_status_pkey PRIMARY KEY (id);


--
-- Name: deployment_directory deployment_directory_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.deployment_directory
    ADD CONSTRAINT deployment_directory_pkey PRIMARY KEY (id);


--
-- Name: project_link project_link_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.project_link
    ADD CONSTRAINT project_link_pkey PRIMARY KEY (id);


--
-- Name: project project_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.project
    ADD CONSTRAINT project_pkey PRIMARY KEY (id, visibility);


--
-- Name: dataset_artifact dataset_artifact_dataset_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.dataset_artifact
    ADD CONSTRAINT dataset_artifact_dataset_id_fkey FOREIGN KEY (dataset_id) REFERENCES public.dataset(id);


--
-- Name: dataset dataset_collection_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.dataset
    ADD CONSTRAINT dataset_collection_id_fkey FOREIGN KEY (collection_id, collection_visibility) REFERENCES public.project(id, visibility);


--
-- Name: dataset_processing_status dataset_processing_status_dataset_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.dataset_processing_status
    ADD CONSTRAINT dataset_processing_status_dataset_id_fkey FOREIGN KEY (dataset_id) REFERENCES public.dataset(id);


--
-- Name: project_link project_link_collection_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.project_link
    ADD CONSTRAINT project_link_collection_id_fkey FOREIGN KEY (collection_id, collection_visibility) REFERENCES public.project(id, visibility);


--
-- PostgreSQL database dump complete
--

