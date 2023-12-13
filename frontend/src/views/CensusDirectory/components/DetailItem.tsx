import Link from "next/link";
import { ItemContainer, ItemLabel } from "../styles";

function formatValueIfNumber(value: string | number) {
  const formatter = new Intl.NumberFormat("en", {
    notation: "compact",
    compactDisplay: "short",
    maximumFractionDigits: 1,
  });
  if (typeof value === "number") {
    return formatter.format(value);
  }
  return value;
}

function DetailItem(props: {
  label: string;
  children?: string | number;
  link?: string;
  onClick?: () => void;
}) {
  return props.children ? (
    <ItemContainer>
      <ItemLabel>{props.label}</ItemLabel>
      {props.link ? (
        <Link href={props.link} passHref>
          <a onClick={props.onClick}> {formatValueIfNumber(props.children)}</a>
        </Link>
      ) : (
        formatValueIfNumber(props.children)
      )}
    </ItemContainer>
  ) : null;
}

export default DetailItem;
