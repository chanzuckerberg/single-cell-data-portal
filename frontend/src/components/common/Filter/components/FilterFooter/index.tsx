import { FOOTER_COMPONENT_MAP } from "./utils";

export const FilterFooter = ({
  componentId,
}: {
  componentId?: keyof typeof FOOTER_COMPONENT_MAP;
}) => {
  if (!componentId) {
    return null;
  }

  const Component = FOOTER_COMPONENT_MAP[componentId] || "Hello";

  return <Component />;
};
