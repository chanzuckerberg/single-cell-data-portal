interface IconProps {
  className: string;
}

const CellxgeneIconSmall = (props: IconProps): JSX.Element => {
  const { className } = props;

  return (
    <svg
      className={className}
      width="20"
      height="20"
      viewBox="0 0 20 20"
      fill="none"
      xmlns="http://www.w3.org/2000/svg"
    >
      <path
        fillRule="evenodd"
        clipRule="evenodd"
        d="M20 0V20H0V0H20ZM16.6667 7.5H7.5V16.6667H16.6667V7.5ZM5.83333 7.5H3.33333V16.6667H5.83333V7.5ZM16.6667 3.33333H7.5V5.83333H16.6667V3.33333Z"
        fill="#0073ff"
      />
      <path d="M14.1667 10H10V14.1667H14.1667V10Z" fill="#0073ff" />
    </svg>
  );
};

export default CellxgeneIconSmall;
