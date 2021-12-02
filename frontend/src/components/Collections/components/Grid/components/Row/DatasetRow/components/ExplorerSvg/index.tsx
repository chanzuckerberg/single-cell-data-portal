function SvgComponent(props: unknown): JSX.Element {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      width={16}
      height={16}
      fill="currentColor"
      {...props}
    >
      <path
        fillRule="evenodd"
        d="M5.143 5.143h9.286v9.286H5.143V5.143zm7.143 2.143h-5v5h5v-5z"
      />
      <path d="M1.571 5.143h2.143v9.286H1.571V5.143zm3.572-3.572h9.286v2.143H5.143V1.571z" />
    </svg>
  );
}

export default SvgComponent;
